'use strict';
/* =====================================================================
   fem_sim.js — FEM LAB "FEM Simulation" katmanı (2B, veri güdümlü yapım animasyonu)

   Neden ayrı dosya: simülasyon bir çözücü DEĞİLDİR. Çözücünün (index.html CORE) ürettiği ağı, a_j katsayılarını,
   u_h(x) değerlerini ve referansı; basis_core.js'in (çözücünün shapeFunctions'ından) ürettiği N_k örneklerini TÜKETİR.
   Burada hiçbir denklem kurulmaz, çözülmez, hiçbir katsayı uydurulmaz; eksik veri 0 ile değil null ile temsil edilir.
   DOM'a dokunmaz: veri sözleşmesi, aşama listesi ve SVG metni üretir (node'da test edilebilsin diye). Düğmeler,
   zamanlayıcı ve durum (state.sim) index.html'dedir; tek durum nesnesi kuralı orada korunur.

   Durağan problemde animasyon ilerlemesi FİZİKSEL ZAMAN DEĞİLDİR, yapım adımıdır; temporal alanı null kalır.
   Zamana bağlı çözücü eklendiğinde temporal = { times, snapshots } gerçek adımlarla doldurulur; o zamana kadar
   physicalTimeAvailable false'tur ve hiçbir yerde "zaman" etiketi kullanılmaz.

   Bölüm haritası (bir bölüm yalnız üstündekileri kullanır):
     1. VERİ SÖZLEŞMESİ   buildSimulationData
     2. AŞAMALAR          buildStages, locateStep, frameState, valueRange (2B ve 3B için ortak)
     3. YARDIMCILAR       fmt, niceTicks, conditionLabel
     4. ÇİZİM             renderFrame
     5. SORGU             probeSample
   ===================================================================== */

/* ---------- Radyan ekseni ----------
   Neden üst düzeyde: ana sayfanın genel grafik çizicisi (drawChart) de aynı işaretleri kullansın; tek tanım. Trigonometrik
   girdide (sin, cos, tan) x bir açıdır ve eksen π katlarıyla yazılır: 0, π/2, π, 3π/2 … Yalnız işaretin YAZILIŞI değişir;
   çizilen değerler ve örnekleme aynıdır (x zaten radyandır, trigonometrik fonksiyonlar radyan alır). */

/* Girdi metinlerinde "açı" anlamında trigonometrik çağrı var mı? sin(…)/cos(…)/tan(…) çağrısının argümanı π içermiyorsa x
   açıdır; sin(pi*x) gibi argümanda x boyutsuz bir orandır (π zaten içeride), eksen π katlarına çevrilmez. sinh/cosh/tanh/asin
   açı değildir ve sayılmaz. İç içe parantez sayılarak okunur. */
function usesTrigonometry(texts) {
  for (const raw of texts) {
    const text = String(raw ?? ''), re = /(?<![A-Za-z0-9_])(sin|cos|tan)\s*\(/gi;
    let m;
    while ((m = re.exec(text))) {
      let depth = 1, i = re.lastIndex;
      while (i < text.length && depth) { depth += text[i] === '(' ? 1 : text[i] === ')' ? -1 : 0; i++; }
      if (!/(?<![A-Za-z0-9_])pi(?![A-Za-z0-9_])/i.test(text.slice(re.lastIndex, i - 1))) return true;
    }
  }
  return false;
}

/* π'nin kesirli katları: π/6, π/4, π/3, π/2, π, 2π, 4π … (π/12 gibi daha ince adım yok: dar aralıkta okunmaz, onluğa dönülür) ≤ count işaret veren EN KÜÇÜK adım seçilir; aralıkta üçten az
   işaret çıkarsa null (çağıran onluk işaretlere döner: örn. [0, 1]). Dönüş: [{ value, label }]. */
function angularTicks(lo, hi, count = 5) {
  if (!(hi > lo) || !Number.isFinite(lo + hi)) return null;
  const STEPS = [[1, 6], [1, 4], [1, 3], [1, 2], [1, 1], [2, 1], [4, 1], [5, 1], [10, 1], [20, 1], [50, 1], [100, 1]];
  const gcd = (a, b) => (b ? gcd(b, a % b) : Math.abs(a));
  for (const [p, q] of STEPS) {
    const step = Math.PI * p / q, eps = 1e-9 * (hi - lo);
    const kLo = Math.ceil((lo - eps) / step), kHi = Math.floor((hi + eps) / step);
    const n = kHi - kLo + 1;
    if (n > count) continue;
    if (n < 3) return null;
    const out = [];
    for (let k = kLo; k <= kHi; k++) {
      const num = k * p, g = gcd(num, q) || 1, a = num / g, d = q / g;
      const label = num === 0 ? '0' : `${a < 0 ? '−' : ''}${Math.abs(a) === 1 ? '' : Math.abs(a)}π${d === 1 ? '' : '/' + d}`;
      out.push({ value: k * step, label });
    }
    return out;
  }
  return null;
}

/* ---------- Çözüm kâğıdı gezintisi ----------
   Neden: Solution paneli kâğıdın adım adım görselleştirmesidir; her adımda kâğıdın o bölümünün KENDİ blokları gösterilir (ikinci
   bir türetme sistemi yok). Kâğıt sunucudan gelen markdown metnidir: "## N. Başlık" bölümleri ve $$…$$ blokları. Bölüm numarası
   dinamik olduğu için bölümler başlık metninden tanınır. Bir bölüm ya da blok bulunamazsa null döner ve çağıran yerleşik formüle
   düşer (uydurma ya da eski bir kâğıdın bloğu asla gösterilmez: çağıran anahtarı kontrol eder). */
function parseSheet(text) {
  const out = {};
  const kind = (t) => (/^Given\b/i.test(t) ? 'given' : /^\d+\.\s*Notation/i.test(t) ? 'notation' : /^\d+\.\s*Mesh\b/i.test(t) ? 'mesh' : /^\d+\.\s*Shape functions/i.test(t) ? 'shape'
    : /Reconstruction/i.test(t) ? 'recon' : /^\d+\.\s*Approximate solution/i.test(t) ? 'approx' : /Galerkin Weighted/i.test(t) ? 'galerkin' : /Weak Formulation/i.test(t) ? 'weak'
    : /Element matrices/i.test(t) ? 'element' : /^\d+\.\s*Assembly/i.test(t) ? 'assembly' : /Boundary conditions/i.test(t) ? 'bc'
    : /^\d+\.\s*Solution\b/i.test(t) ? 'solution' : /Error measures/i.test(t) ? 'error' : null);   // 'Reconstruction of the approximate solution' önce: başlık 'Approximate solution' içerir
  for (const part of String(text || '').split(/^## /m).slice(1)) {
    const title = part.split('\n', 1)[0].trim(), key = kind(title);
    if (key) out[key] = { title, blocks: [...part.matchAll(/\$\$([\s\S]*?)\$\$/g)].map((m) => m[1].trim()) };
  }
  return Object.keys(out).length ? out : null;
}

/* Bir adımın kâğıt blokları: { tex: [...], source } ya da null. ctx: { order, p, n }. Eleman numarasına bağlı bloklar (K⁽ᵉ⁾_ij, F⁽ᵉ⁾_i,
   u_h⁽ᵉ⁾, R_e) için kâğıdın gösterdiği elemanlardan seçili elemana EN YAKIN olan alınır: sabit katsayıda kâğıt yalnız e = 1'i
   açar ("tüm elemanlar için aynı"), n > 8'de ilk iki ve son elemanı. */
function sheetPick(sections, stage, fs, ctx) {
  if (!sections) return null;
  const sec = (key) => sections[key] || null, blocks = (key) => (sec(key) ? sec(key).blocks : []);
  const { selEl } = fs, { order, p, n } = ctx, k = fs.local;
  const nearest = (list) => (list.length ? list.reduce((b, e) => (Math.abs(e - selEl) < Math.abs(b - selEl) ? e : b)) : null);
  const esc = (t) => t.replace(/[.*+?^${}()|[\]\\]/g, '\\$&');
  /* make(e) → RegExp; seçili elemana en yakın, bloğu olan elemanın bloğu */
  const forElement = (bl, make, prefer = selEl) => {
    const have = [];
    for (let e = 1; e <= n; e++) if (bl.some((b) => make(e).test(b))) have.push(e);
    if (!have.length) return null;
    const e = have.reduce((b, q) => (Math.abs(q - prefer) < Math.abs(b - prefer) ? q : b));
    return { e, block: bl.find((b) => make(e).test(b)) };
  };
  const has = (bl, sub) => bl.filter((b) => b.includes(sub));
  const starts = (bl, sub) => bl.filter((b) => b.startsWith(sub));
  const done = (tex, title) => { const t = tex.filter(Boolean); return t.length ? { tex: t, source: title } : null; };
  /* Panel 1 (Approximation): kâğıdın Given / Notation / Mesh / Shape functions / Approximate solution bölümleri. Kâğıdın referans elemanı
     s = x − x_e ∈ [0, h]'dir (Jacobian 1); ekrandaki çizimler de aynı s ile etiketlenir. Aşama kimlikleri çözüm paneliyle çakışmaz. */
  const afterLast = (bl, i) => (i >= 0 && i < bl.length ? bl[i] : null);
  switch (stage.id) {
    case 'domain': {
      const bl = blocks('given'); if (!bl.length) return null;
      const iEq = bl.findIndex((b) => b.includes('\\in \\Omega = (')), L = afterLast(bl, iEq + 1), cond = afterLast(bl, iEq + 2);
      if (iEq < 0 || !L || !L.startsWith('\\mathcal{L}(u)') || !cond) return null;
      return done(k === 0 ? [L] : k === 1 ? [bl[iEq]] : [cond], sec('given').title);
    }
    case 'mesh': {
      const g = blocks('given'), bl = blocks('mesh'); if (!bl.length) return null;
      const params = starts(g, 'n = ').slice(0, 1), last = fs.finalMesh;
      return done(k === 0 ? [...params, ...starts(g, 'N = np + 1').slice(0, 1)] : last ? [bl[0], bl[1]] : [bl[0]], k === 0 ? sec('given').title : sec('mesh').title);
    }
    case 'element': {
      const sh = blocks('shape'), bl = blocks('mesh'); if (!sh.length) return null;
      return done(k === 0 ? [sh[0]] : [...starts(sh, 'x = x_e + s').slice(0, 1), bl[1]], sec('shape').title);
    }
    case 'local': {
      const sh = blocks('shape'); if (!sh.length) return null;
      const all = starts(sh, '\\begin{aligned} N_{1}(s)')[0] || null, rows = all ? all.replace(/^\\begin\{aligned\}|\\end\{aligned\}$/g, '').split('\\\\').map((r) => r.trim()).filter(Boolean) : [];
      const own = rows.length === p + 1 ? `\\begin{aligned} ${rows[k]} \\end{aligned}` : all;
      return done([own, ...has(sh, '\\big[N_i(s_j)\\big]').slice(0, 1)], sec('shape').title);
    }
    case 'global': {
      const ap = blocks('approx'), no = blocks('notation'), sh = blocks('shape'); if (!ap.length && !no.length) return null;
      if (k === 0) return done([ap[0], ap[1]], sec('approx').title);
      if (k === 1) return done(starts(no, 'N_j(x) = \\begin{cases}').slice(0, 1), sec('notation').title);
      return done(has(sh, 'continuity:').slice(0, 1), sec('shape').title);
    }
    case 'ansatz': {
      const bl = blocks('approx'); if (!bl.length) return null;
      const el = forElement(bl, (e) => new RegExp(`^u_h\\^\\{\\(${e}\\)\\}\\(x\\)`));
      return done(k === 0 ? [bl[0], bl[1]] : [bl[2], el && el.block], sec('approx').title);
    }
    case 'residual': {
      const bl = blocks('galerkin'); if (!bl.length) return null;
      if (k === 0) { const el = forElement(bl, (e) => new RegExp(`^R_\\{${e}\\}\\(x\\)`)); return done([bl[0], el && el.block], sec('galerkin').title); }
      return done([...has(bl, '\\phi_i(x)\\,R(u_h'), ...has(bl, '\\phi_i(x) = N_i(x)'), ...has(bl, 'removed:').slice(0, 2)], sec('galerkin').title);
    }
    case 'weak': {
      const bl = blocks('weak'); if (!bl.length) return null;
      return done(k === 0 ? [bl[0], bl[1]] : [bl[2], ...has(bl, '\\underbrace').slice(0, 1), bl[bl.length - 1]], sec('weak').title);
    }
    case 'elementK': {
      const bl = blocks('element'); if (!bl.length) return null;
      const title = sec('element').title, kind = fs.kPhase;
      if (kind === 'intro') return done([...starts(bl, 'K^{(e)}_{ij} = \\int_{x_e}').slice(0, 1), ...starts(bl, 'x = x_e + s').slice(0, 1), ...starts(bl, 'K^{(e)}_{ij} = \\int_0^h').slice(0, 1)], title);
      if (kind === 'K') {
        const { i, j } = fs.cell, cellRe = (e) => new RegExp(`^\\\\begin\\{aligned\\} K\\^\\{\\(${e}\\)\\}_\\{${i + 1}${j + 1}\\}`);
        const c = forElement(bl, cellRe), head = c && bl.find((b) => new RegExp(`^e = ${c.e}:`).test(b));
        return done([head, c && c.block], title);
      }
      if (kind === 'F') {
        const { fIndex } = fs, c = forElement(bl, (e) => new RegExp(`^\\\\begin\\{aligned\\} F\\^\\{\\(${e}\\)\\}_\\{l,${fIndex + 1}\\}`));
        return done([c && c.block], title);
      }
      const km = forElement(bl, (e) => new RegExp(`^\\\\mathbf\\{K\\}\\^\\{\\(${e}\\)\\}`)), fm = forElement(bl, (e) => new RegExp(`^\\\\mathbf\\{F\\}_l\\^\\{\\(${e}\\)\\}`));
      return done([km && km.block, fm && fm.block], title);
    }
    case 'assembly': {
      const bl = blocks('assembly'); if (!bl.length) return null;
      const out = [bl[0]];
      if (k >= 2) out.push(bl[1], bl[2]);
      if (k === n) out.push(...starts(bl, '\\mathbf{K} = ').slice(0, 1));
      return done(out, sec('assembly').title);
    }
    case 'bc': {
      const a8 = blocks('assembly'), a9 = blocks('bc'); if (!a9.length) return null;
      const matrices = starts(a8, '\\mathbf{K} = ').slice(0, 1), fb = starts(a8, '\\mathbf{F}_b').slice(0, 1);
      const essential = a9.slice(0, 3), reduced = [a9[3], a9[a9.length - 1]].filter((b, i, arr) => arr.indexOf(b) === i);
      const parts = order === 2 ? [matrices, fb, essential, reduced] : [matrices, essential, reduced];
      return done(parts[Math.min(k, parts.length - 1)], `${a8.length ? sec('assembly').title + ' → ' : ''}${sec('bc').title}`);
    }
    case 'solve': {
      const a9 = blocks('bc'), a10 = blocks('solution'); if (!a9.length && !a10.length) return null;
      return done(k === 0 ? [a9[a9.length - 1]] : a10.slice(0, 3), k === 0 ? sec('bc').title : sec('solution').title);
    }
    case 'approx': {
      const bl = blocks('recon'); if (!bl.length) return null;
      const el = forElement(bl, (e) => new RegExp(`^u_h\\^\\{\\(${e}\\)\\}\\(x\\)`));
      return done([bl[0], el && el.block], sec('recon').title);
    }
    case 'reference': {
      const bl = blocks('error'); if (!bl.length || !bl[0].startsWith('u(x)')) return null;
      return done([bl[0]], sec('error').title);
    }
    case 'error': {
      const bl = blocks('error'); if (!bl.length) return null;
      return done(k === 0 ? [bl[1]] : [...starts(bl, '\\max_j'), ...starts(bl, '\\max_x'), ...starts(bl, '\\|u - u_h\\|_{L^2} = '), ...starts(bl, '|u - u_h|_{H^1} = ')].slice(0, 4), sec('error').title);
    }
    default: return null;
  }
}

/* Neden fabrika: basis_core.js gibi, çekirdek fonksiyonları dışarıdan alır (index.html ile node testleri aynı kodu
   çalıştırır). deps.basis = femExtension(...) sonucu; deps.evaluate = CORE evaluate(mesh, c, x). */
function femSimulation(deps) {
  const { basis, evaluate, elementSystem } = deps;   // elementSystem: çekirdeğin eleman matrisi fonksiyonu (yoksa cebir aşamaları yok)
  if (!basis || typeof evaluate !== 'function') throw new Error('femSimulation needs basis_core and the solver\'s evaluate.');
  const SAMPLES = 33;   // eleman başına örnek (uçlar dahil): P4'ün en fazla 4 salınımı da düzgün çizilsin

  /* ---------- 1. VERİ SÖZLEŞMESİ ----------
     Neden tek nesne: 2B çizici (bugün) ve 3B çizici (Stage D) aynı veriyi okusun; çizici hesap yapmasın.
     Zorunlu: problem, discretization, basis, solution.nodalCoefficients/samples.fem. İsteğe bağlı (yoksa null):
     solution.samples.reference/error, metadata.maxError. Desteklenmeyen (bugün): temporal (null). */
  /* Cebir verisi (Solution paneli için). Neden burada hesap yok: eleman matrisleri ÇEKİRDEĞİN elementSystem'idir; global K, F ve
     a çözücünün sonucudur. Burada yalnız defter tutma var: (1) birinci mertebede çekirdeğin yaptığı (−1) işaret çevirmesi,
     (2) sınır akısı F_b, (3) hangi düğümlerin sabit, hangi satır/sütunların kullanıldığı (çekirdekteki kuralın aynısı).
     tests_sim.js bunların çekirdekle tuttuğunu kanıtlar: ΣK⁽ᵉ⁾ = K, ΣF⁽ᵉ⁾ + F_b = F, indirgenmiş sistemin çözümü = a. */
  function buildAlgebra(problem, mesh, c, solved, local, elementX) {
    if (typeof elementSystem !== 'function' || !solved || !Array.isArray(solved.K) || !Array.isArray(solved.F)) return null;
    const { n, p, h, dofs } = mesh, sign = problem.order === 1 ? -1 : 1, J = 2 / h, S = local.xi.length;
    const flip = (v) => (v === 0 ? 0 : sign * v);
    const elementK = [], elementF = [];
    for (let e = 0; e < n; e++) {
      const { Ke, Fe } = elementSystem(problem, mesh, e, 'fine');
      elementK.push(Ke.map((row) => row.map(flip)));
      elementF.push(Fe.map(flip));
    }
    /* Integrand terimleri (grafik için): N_i' α N_j', −N_i β N_j', −N_i γ N_j. Alan = K⁽ᵉ⁾_ij (çekirdek kuadratürüyle). */
    const terms = new Float64Array(n * (p + 1) * (p + 1) * 3 * S);
    for (let e = 0; e < n; e++) for (let q = 0; q < S; q++) {
      const x = elementX[e][q], av = problem.coef.alpha.f(x), bv = problem.coef.beta.f(x), gv = problem.coef.gamma.f(x);
      for (let i = 0; i <= p; i++) for (let j = 0; j <= p; j++) {
        const base = (((e * (p + 1) + i) * (p + 1) + j) * 3) * S + q;
        terms[base] = flip(local.dN[i][q] * J * av * local.dN[j][q] * J);
        terms[base + S] = flip(-local.N[i][q] * bv * local.dN[j][q] * J);
        terms[base + 2 * S] = flip(-local.N[i][q] * gv * local.N[j][q]);
      }
    }
    const fterms = new Float64Array(n * (p + 1) * S);          // yük integrandı: −N_i f (çekirdeğin işaretiyle)
    for (let e = 0; e < n; e++) for (let q = 0; q < S; q++) {
      const fv = problem.coef.f.f(elementX[e][q]);
      for (let i = 0; i <= p; i++) fterms[(e * (p + 1) + i) * S + q] = flip(-local.N[i][q] * fv);
    }
    const Fl = new Array(dofs).fill(0), Fb = new Array(dofs).fill(0);
    for (let e = 0; e < n; e++) for (let i = 0; i <= p; i++) Fl[e * p + i] += elementF[e][i];
    if (problem.order === 2) {
      const flux = (side, x) => {
        const { type, value, slope } = problem.bc[side], alpha = problem.coef.alpha.f(x);
        return type === 'neumann' ? alpha * value : type === 'flux' ? value : type === 'both' ? alpha * slope : 0;
      };
      Fb[0] -= flux('left', mesh.a);
      Fb[dofs - 1] += flux('right', mesh.b);
    }
    const fixed = [];
    if (['dirichlet', 'both'].includes(problem.bc.left.type)) fixed.push([0, problem.bc.left.value]);
    if (['dirichlet', 'both'].includes(problem.bc.right.type)) fixed.push([dofs - 1, problem.bc.right.value]);
    const fixedSet = new Set(fixed.map(([i]) => i)), free = [];
    for (let i = 0; i < dofs; i++) if (!fixedSet.has(i)) free.push(i);
    const rows = problem.bc.left.type === 'both' ? Array.from({ length: dofs - 1 }, (_, i) => i)
      : problem.bc.right.type === 'both' ? Array.from({ length: dofs - 1 }, (_, i) => i + 1) : free;
    const K = solved.K, F = solved.F;
    const reducedK = rows.map((i) => free.map((j) => K[i][j]));
    const reducedF = rows.map((i) => { let v = F[i]; for (const [j, value] of fixed) v -= K[i][j] * value; return v; });
    return { elementK, elementF, terms, fterms, termSamples: S, Fl, Fb, K, F, fixed, free, rows, reducedK, reducedF,
      a: c.slice(), condition: Number.isFinite(solved.condition) ? solved.condition : null };
  }
  const termAt = (A, p, e, i, j, t, s) => A.terms[((((e * (p + 1) + i) * (p + 1) + j) * 3 + t) * A.termSamples) + s];

  function buildSimulationData({ problem, mesh, c, xs, reference = null, maxError = null, angular = false, solved = null, tex = null }) {
    if (!mesh || !(mesh.n >= 1) || !(mesh.p >= 1)) throw new Error('A mesh is required.');
    if (!Array.isArray(c) || c.length !== mesh.dofs || !c.every(Number.isFinite)) throw new Error(`${mesh.dofs} finite nodal coefficients are required.`);
    if (!Array.isArray(xs) || xs.length < 2) throw new Error('Sample points are required.');
    const { n, p, h, a, b, dofs } = mesh;
    const coordinates = Array.from({ length: dofs }, (_, i) => mesh.nodeX(i));          // çözücünün kendi düğüm koordinatı
    const connectivity = Array.from({ length: n }, (_, e) => Array.from({ length: p + 1 }, (_, k) => e * p + k + 1));   // 1 tabanlı global düğüm
    const elementBoundaries = Array.from({ length: n + 1 }, (_, e) => (e === n ? b : a + e * h));
    const local = basis.referenceSamples(p, SAMPLES);                                    // N_k(ξ): çözücünün shapeFunctions'ı
    const elementX = connectivity.map((_, e) => local.xi.map((t) => a + e * h + h * (1 + t) / 2));
    const support = Array.from({ length: dofs }, () => []);                              // N_j'nin sıfırdan farklı olduğu elemanlar
    connectivity.forEach((nodes, e) => nodes.forEach((I, k) => support[I - 1].push({ element: e + 1, local: k + 1 })));
    const algebra = buildAlgebra(problem, mesh, c, solved, local, elementX);
    const fem = xs.map((x) => evaluate(mesh, c, x));
    const finiteOrNull = (v) => (Number.isFinite(v) ? v : null);
    const refValues = reference ? xs.map((x) => finiteOrNull(reference.u(x))) : null;
    const errValues = refValues ? refValues.map((r, i) => (r === null || !Number.isFinite(fem[i]) ? null : r - fem[i])) : null;
    /* Girdi fonksiyonları: çözücünün denklemi kurarken kullandığı AYNI ayrıştırılmış ifadeler (problem.coef) örnek
       noktalarında değerlendirilir. Sıfır olan α, β, γ terimi denklemde yoktur, çizilmez; f her zaman çizilir (f = 0 da bilgidir). */
    const INPUT_NAMES = { alpha: 'α(x)', beta: 'β(x)', gamma: 'γ(x)', f: 'f(x)' };
    const inputFunctions = ['alpha', 'beta', 'gamma', 'f']
      .filter((key) => problem.coef?.[key] && (key === 'f' || problem.constant?.[key] !== 0))
      .map((key) => ({ key, name: INPUT_NAMES[key], constant: Number.isFinite(problem.constant?.[key]) ? problem.constant[key] : null,
        values: xs.map((x) => finiteOrNull(problem.coef[key].f(x))) }));
    const bc = (side, x) => {
      const cnd = problem.bc?.[side] || { type: 'none' };
      return { side, x, type: cnd.type, value: Number.isFinite(cnd.value) ? cnd.value : null, slope: Number.isFinite(cnd.slope) ? cnd.slope : null };
    };
    return {
      version: 1,
      problem: {
        equationClass: problem.order === 2 ? 'stationary_second_order' : 'stationary_first_order',
        equationOrderX: problem.order, equationOrderT: 0, domain: [a, b], timeDomain: null,
        initialValueProblem: !!problem.ivp, boundaryConditions: [bc('left', a), bc('right', b)], initialConditions: null, inputFunctions,
        angularX: !!angular,  // trigonometrik girdi: x ekseni radyan (π katları) yazılır
        tex: tex || null      // denklem ve koşulların LaTeX'i (index.html'deki çekirdek yazıcısından), formül altyazıları için
      },
      discretization: { elementCount: n, nodeCount: dofs, polynomialDegree: p, h, coordinates, connectivity, elementBoundaries },
      basis: { referenceNodes: local.nodes, xi: local.xi, N: local.N, dN: local.dN, elementX, support },
      algebra,
      solution: { nodalCoefficients: c.slice(), samples: { x: xs.slice(), fem, reference: refValues, error: errValues } },
      temporal: null,
      metadata: {
        referenceType: reference ? reference.type : null,                 // 'analytical' | 'fine_mesh' | null
        referenceLevel: reference ? reference.level ?? null : null,
        errorAvailable: !!reference,
        errorKind: reference ? (reference.type === 'analytical' ? 'exact' : 'difference') : null,
        physicalTimeAvailable: false,
        maxError: maxError && Number.isFinite(maxError.value) ? { x: maxError.x, value: maxError.value } : null
      }
    };
  }

  /* ---------- 2. AŞAMALAR ----------
     Neden veriden kurulur: aşama ve adım sayısı ağa bağlıdır (n eleman, n·p + 1 terim); referans yoksa karşılaştırma ve
     hata aşamaları hiç eklenmez (boş grafik ya da uydurma hata gösterilmesin). Her adım AYRIK bir durumdur: ara değer
     enterpolasyonu yok, yalnız gerçek nesnelerin (eleman, düğüm, terim) sırayla görünmesi. */
  function buildStages(data) {
    const { elementCount: n, nodeCount, polynomialDegree: p } = data.discretization;
    const many = (k, word) => `${k} ${word}${k === 1 ? '' : 's'}`;
    const order = data.problem.equationOrderX, A = data.algebra;
    /* İki panel: 'approx' = nasıl yaklaşıyoruz (ağ ve taban; denklemden bağımsız), 'solution' = katsayılar nasıl bulunuyor
       (Galerkin cebiri; denkleme özgü) ve sonuç. Aşama sırası FEM'in gerçek sırasıdır; her panelin kendi zaman çubuğu var. */
    const stages = [
      { panel: 'approx', id: 'domain', title: 'Continuous problem', steps: 3,
        text: 'The differential problem as entered: the input functions α, β, γ, f, then the interval and its boundary conditions; nothing is discretized yet.' },
      { panel: 'approx', id: 'mesh', title: 'Mesh generation', steps: Math.min(n, 16) + 1,
        text: `The interval is divided into ${many(n, 'element')} of equal size h; the ${nodeCount} nodes carry the unknowns.` },
      { panel: 'approx', id: 'element', title: 'Element structure', steps: 2,
        text: `Each element is the reference element s ∈ [0, h] shifted to x = x_e + s, and has p + 1 = ${p + 1} local nodes.` },
      { panel: 'approx', id: 'local', title: 'Local shape functions', steps: p + 1,
        text: `On the selected element the ${p + 1} Lagrange shape functions N_k equal 1 at their own node and 0 at the other local nodes.` },
      { panel: 'approx', id: 'global', title: 'Global basis construction', steps: 3,
        text: 'Local shape functions that share a node join into one global basis function N_j, nonzero only on its support.' }
    ];
    /* Solution paneli = çözüm kâğıdının bölümleri sırayla (başlıklar kâğıttakiyle aynı). Metinler tek cümle: asıl anlatım
       her adımın altındaki kâğıt blokları. */
    if (A) {
      stages.push(
        { panel: 'solution', id: 'ansatz', title: 'Approximate solution', steps: 2,
          text: 'u_h = Σ a_j N_j: one unknown coefficient a_j per node.' },
        { panel: 'solution', id: 'residual', title: 'Galerkin weighted residual', steps: 2,
          text: 'The residual R = 𝓛(u_h) − f is made orthogonal to the weight functions w_i = N_i.' });
      if (order === 2) stages.push({ panel: 'solution', id: 'weak', title: 'Weak formulation', steps: 2,
        text: 'Integration by parts moves one derivative onto the weight function.' });
      stages.push(
        { panel: 'solution', id: 'elementK', title: 'Element matrices and load vectors', steps: (p + 1) * (p + 1) + (p + 1) + 2,
          text: 'Each entry is one integral over the element: K⁽ᵉ⁾_ij pairs weight N_i with unknown N_j; F⁽ᵉ⁾_i uses the source f.' },
        { panel: 'solution', id: 'assembly', title: 'Assembly', steps: n + 1,
          text: 'Each element matrix is added at the rows and columns of its global node numbers.' },
        { panel: 'solution', id: 'bc', title: 'Boundary conditions and reduced system', steps: order === 2 ? 4 : 3,
          text: 'Known coefficients move to the right-hand side; their rows and columns are removed.' },
        { panel: 'solution', id: 'solve', title: 'Solution', steps: 2,
          text: 'The reduced system K_red a = F_red is solved for the unknown coefficients.' });
    }
    stages.push(
      { panel: 'solution', id: 'approx', title: 'Reconstruction of the approximate solution', steps: Math.min(nodeCount, 24) + 1,
        text: 'The solved a_j scale the basis functions; their sum is u_h(x).' });
    if (data.metadata.referenceType) {
      const analytical = data.metadata.referenceType === 'analytical';
      stages.push({ panel: 'solution', id: 'reference', title: analytical ? 'Analytical solution comparison' : 'Reference comparison', steps: 2,
        text: analytical ? 'The analytical solution is drawn with the FEM approximation.'
          : 'No analytical solution: the fine-mesh reference (n = 128, p = 4) is drawn with the FEM approximation.' });
      stages.push({ panel: 'solution', id: 'error', title: analytical ? 'Error measures' : 'Difference from the reference', steps: 2,
        text: analytical ? 'e(x) = u(x) − u_h(x): positive where u_h lies below the analytical solution.'
          : 'd(x) = u_ref(x) − u_h(x): not the exact error, which needs an analytical solution.' });
    }
    let offset = 0;
    for (const s of stages) { s.offset = offset; offset += s.steps; }
    return stages;
  }
  /* Bir panelin aşamaları, kendi içinde yeniden numaralanmış (offset 0'dan): her panelin zaman çubuğu bağımsızdır. */
  function stagesFor(stages, panel) {
    let offset = 0;
    return stages.filter((st) => st.panel === panel).map((st) => { const copy = { ...st, offset }; offset += st.steps; return copy; });
  }

  /* Genel adım → (aşama, aşama içi adım). Aralık dışı adımlar kırpılır: durum her zaman geçerli bir kareye karşılık gelsin. */
  function locateStep(stages, step) {
    const total = stages.reduce((s, st) => s + st.steps, 0);
    const k = Math.max(0, Math.min(total - 1, Math.floor(Number.isFinite(step) ? step : 0)));
    const index = stages.findIndex((st) => k < st.offset + st.steps);
    return { index, stage: stages[index], local: k - stages[index].offset, step: k, total };
  }

  /* Karenin anlamsal durumu: hangi aşama, kaç eleman/terim görünür, seçili düğüm/eleman, komşular. Neden ayrı: 2B (SVG) ve
     3B (Three.js) çizici AYNI sayıları kullansın; iki çizici farklı şey göstermesin (tests_sim.js ikisini karşılaştırır). */
  function frameState(data, stages, step, opts = {}) {
    const where = locateStep(stages, step), { stage, local: k } = where;
    const { elementCount: n, nodeCount, polynomialDegree: p, connectivity: C } = data.discretization;
    const selNode = Number.isInteger(opts.node) && opts.node >= 1 && opts.node <= nodeCount ? opts.node : (n >= 2 ? Math.floor(n / 2) * p + 1 : 1);
    const selEl = Number.isInteger(opts.element) && opts.element >= 1 && opts.element <= n ? opts.element : Math.ceil(n / 2);
    const last = k === stage.steps - 1;
    const shownElements = stage.id === 'mesh' && !last ? Math.ceil((k + 1) * n / (stage.steps - 1)) : n;
    const terms = stage.id === 'approx' ? (last ? nodeCount : Math.ceil((k + 1) * nodeCount / (stage.steps - 1))) : 0;
    const newestFrom = stage.id === 'approx' && !last ? Math.ceil(k * nodeCount / (stage.steps - 1)) + 1 : null;   // bu adımda eklenen ilk terim
    const prevShown = stage.id === 'mesh' && k > 0 && !last ? Math.ceil(k * n / (stage.steps - 1)) : stage.id === 'mesh' && last ? n : 0;
    const prevTerms = stage.id === 'approx' && newestFrom !== null ? newestFrom - 1 : null;          // bir önceki kısmi toplamın terim sayısı
    /* Eleman matrisi aşaması: adım 0 giriş, sonra K⁽ᵉ⁾ hücreleri satır satır, sonra F⁽ᵉ⁾ girdileri, son adım tamamı. Seçili hücre adımdan gelir. */
    let kPhase = null, cell = null, fIndex = null;
    if (stage.id === 'elementK') {
      const q = (p + 1) * (p + 1);
      if (k === 0) kPhase = 'intro';
      else if (k <= q) { kPhase = 'K'; cell = { i: Math.floor((k - 1) / (p + 1)), j: (k - 1) % (p + 1) }; }
      else if (k <= q + p + 1) { kPhase = 'F'; fIndex = k - q - 1; }
      else kPhase = 'final';
    }
    const neighbors = new Set();
    for (const { element } of data.basis.support[selNode - 1]) for (const j of C[element - 1]) if (j !== selNode) neighbors.add(j);
    return { ...where, selNode, selEl, shownElements, finalMesh: stage.id !== 'mesh' || last, terms, complete: stage.id === 'approx' && last,
      newestFrom, prevShown, prevTerms, kPhase, cell, fIndex, neighbors: [...neighbors].sort((a, b) => a - b) };
  }
  /* Fonksiyon aşamalarının dikey aralığı: aşama boyunca SABİT (adım ilerlerken eksen zıplamasın), 0 her zaman içeride,
     %8 pay. Değerler yalnız hesaplanmış verilerdir; null/NaN atlanır. Fonksiyon çizmeyen aşamada null. */
  function valueRange(data, stageId) {
    const S = data.solution, B = data.basis, a = S.nodalCoefficients;
    let values;
    if (stageId === 'local' || stageId === 'global') values = B.N.flat().concat([0, 1]);
    else if (stageId === 'approx') {
      values = S.samples.fem.slice();
      B.support.forEach((sup, j) => { for (const { local } of sup) for (const v of B.N[local - 1]) values.push(a[j] * v); });
    } else if (stageId === 'reference') values = S.samples.fem.concat(S.samples.reference || []);
    else if (stageId === 'error') values = S.samples.error || [];
    else return null;
    let lo = Infinity, hi = -Infinity;
    for (const v of values) if (v !== null && Number.isFinite(v)) { lo = Math.min(lo, v); hi = Math.max(hi, v); }
    if (!Number.isFinite(lo)) { lo = -1; hi = 1; }
    lo = Math.min(lo, 0); hi = Math.max(hi, 0);
    if (hi - lo < 1e-300) { lo -= 1; hi += 1; }
    const m = (hi - lo) * 0.08;
    return { lo: lo - m, hi: hi + m };
  }

  /* x ekseni işaretleri: radyan modunda π katları (yeterli işaret çıkarsa), değilse onluk. Dönüş: [{ value, label }]. */
  function xTicks(lo, hi, count, angular) {
    if (angular) { const t = angularTicks(lo, hi, count); if (t) return t; }
    return niceTicks(lo, hi, count).map((value) => ({ value, label: fmt(value) }));
  }

  /* ---------- Formül altyazıları ----------
     Her adımın altında o adımın formülü, KENDİ denklemine göre yazılır (yalnız sıfır olmayan terimler). Sabit kalıplar
     (varsayım, montaj, esas koşul) çözüm kâğıdının bölümleriyle AYNI dizgedir; tests_sim.js kâğıt çıktısında bulunduğunu
     doğrular. Kâğıt tek kaynak: burada ikinci bir türetme sistemi yoktur, yalnız kâğıdın ilgili satırı adım adım gösterilir. */
  function stageFormulas(data, fs, opts = {}) {
    const { stage, local: k, selNode, selEl } = fs;
    const D = data.discretization, P = data.problem, A = data.algebra, order = P.equationOrderX;
    const { elementCount: n, nodeCount: N, polynomialDegree: p } = D;
    const keys = new Set((P.inputFunctions || []).map((f) => f.key));
    const L = [keys.has('alpha') && "(\\alpha u')'", keys.has('beta') && "\\beta u'", keys.has('gamma') && '\\gamma u'].filter(Boolean).join(' + ') || '0';
    const kernel = order === 2
      ? ["N_i'\\,\\alpha\\,N_j'", keys.has('beta') && "- N_i\\,\\beta\\,N_j'", keys.has('gamma') && '- N_i\\,\\gamma\\,N_j'].filter(Boolean).join(' ')
      : ['N_i(' + [keys.has('beta') && "\\beta N_j'", keys.has('gamma') && '\\gamma N_j'].filter(Boolean).join(' + ') + ')'].join('');
    const Fi = order === 2 ? 'F_i = -\\int_\\Omega N_i\\,f\\,dx + F_{b,i}' : 'F_i = \\int_\\Omega N_i\\,f\\,dx';
    const out = [];
    switch (stage.id) {
      case 'domain':
        if (P.tex) out.push(P.tex.equation, P.tex.conditions);
        break;
      case 'mesh': out.push(`h = \\frac{x_1 - x_0}{n} = ${fmt(D.h)},\\qquad x_i = x_0 + \\frac{i\\,h}{p}\\qquad n_n = n\\,p + 1 = ${N}`); break;
      case 'element': out.push(`s = x - x_e \\in [0,\\,h],\\qquad x = x_e + s,\\qquad dx = ds\\qquad (h = ${fmt(D.h)})`, 'j = (e-1)\\,p + k\\qquad \\text{(local node } k \\to \\text{ global node } j)'); break;
      case 'local': out.push('N_k(s) = \\prod_{m \\ne k} \\frac{s - s_m}{s_k - s_m},\\qquad s_k = (k - 1)\\,\\frac{h}{p},\\qquad N_k(s_m) = \\delta_{km}'); break;
      case 'global': out.push('N_j(x)\\big|_{\\Omega_e} = N_k^{(e)}(x),\\quad j = (e-1)\\,p + k,\\qquad N_j = 0 \\text{ elsewhere}'); break;
      case 'ansatz':
        out.push(`u_h(x) = \\sum_{j=1}^{${N}} a_j\\,N_j(x)\\qquad \\text{(} a_j \\text{: nodal coefficients to be determined)}`, 'N_j(x_i) = \\delta_{ij} \\;\\Rightarrow\\; u_h(x_i) = a_i');
        break;
      case 'residual':
        if (k === 0) out.push(`\\mathcal{L}[u] = ${L},\\qquad R(x) = \\mathcal{L}[u_h](x) - f(x) \\ne 0`);
        else out.push(`\\int_\\Omega w_i\\,R\\,dx = 0,\\qquad w_i = N_i,\\qquad i = 1,\\dots,${N}`);
        break;
      case 'weak':
        out.push(k === 0 ? `\\int_\\Omega w_i\\,(\\alpha u_h')'\\,dx = \\Big[w_i\\,\\alpha u_h'\\Big] - \\int_\\Omega w_i'\\,\\alpha u_h'\\,dx`
          : `K_{ij} = \\int_\\Omega ${kernel}\\,dx,\\qquad \\mathbf{K}\\,\\mathbf{a} = \\mathbf{F}`, ...(k === 1 ? [Fi] : []));
        break;
      case 'elementK': {
        const defK = `K^{(e)}_{ij} = \\int_{x_e}^{x_{e+1}} ${kernel}\\,dx`, defF = `F^{(e)}_i = ${order === 2 ? '-' : ''}\\int_{x_e}^{x_{e+1}} N_i\\,f\\,dx`;
        if (fs.kPhase === 'intro' || !A) out.push(defK, defF);
        else if (fs.kPhase === 'K') out.push(defK, `K^{(${selEl})}_{${fs.cell.i + 1}${fs.cell.j + 1}} = ${fmt(A.elementK[selEl - 1][fs.cell.i][fs.cell.j])}`);
        else if (fs.kPhase === 'F') out.push(defF, `F^{(${selEl})}_{${fs.fIndex + 1}} = ${fmt(A.elementF[selEl - 1][fs.fIndex])}`);
        else out.push(`\\mathbf{K}^{(${selEl})},\\ \\mathbf{F}^{(${selEl})}\\ \\text{complete}`);
        break;
      }
      case 'assembly':
        out.push('K_{IJ} = \\sum_{e} K^{(e)}_{ij},\\qquad F_I = \\sum_e F^{(e)}_i,\\qquad I = (e-1)p + i,\\ J = (e-1)p + j');
        break;
      case 'bc': {
        const essential = (A ? A.fixed : []).map(([idx, v]) => `a_{${idx + 1}} = ${fmt(v)}\\ \\text{(essential)}`);
        const hasFlux = order === 2 && k >= 1;
        if (hasFlux) out.push('\\mathbf{F} = \\mathbf{F}_l + \\mathbf{F}_b,\\qquad \\mathbf{F}_b = [-Q_0,\\,0,\\dots,0,\\,Q_1]^T,\\qquad Q = \\alpha\\,u\'');
        if (k >= (order === 2 ? 2 : 1) && essential.length) out.push(essential.join(',\\qquad '));
        if (k >= (order === 2 ? 3 : 2)) out.push('K_{\\mathrm{red}}\\,a_{\\mathrm{free}} = F - K_{\\mathrm{fixed}}\\,a_{\\mathrm{fixed}}');
        break;
      }
      case 'solve':
        out.push(`a_{\\mathrm{free}} = K_{\\mathrm{red}}^{-1}\\,F_{\\mathrm{red}}` + (A && A.condition ? `,\\qquad \\kappa_1 \\approx ${A.condition.toExponential(1)}` : ''));
        if (k >= 1) out.push('a_j = u_h(x_j)\\qquad \\text{(nodal values)}');
        break;
      case 'approx': out.push(fs.complete ? 'u_h(x) = \\sum_{j} a_j\\,N_j(x)' : `\\sum_{j \\le ${fs.terms}} a_j\\,N_j(x)\\qquad \\text{(partial sum, not yet } u_h)`); break;
      case 'reference': out.push('u(x)\\ \\text{and}\\ u_h(x)'); break;
      case 'error': out.push(data.metadata.errorKind === 'exact' ? 'e(x) = u(x) - u_h(x)' : 'd(x) = u_{\\mathrm{ref}}(x) - u_h(x)'); break;
      default: break;
    }
    return out;
  }

  /* ---------- 3. YARDIMCILAR ---------- */
  /* Ekran sayısı: 4 anlamlı basamak; çok küçük/büyükte üstel. Sonsuz/NaN asla sayı gibi yazılmasın. */
  function fmt(v) {
    if (!Number.isFinite(v)) return '—';
    const m = Math.abs(v);
    if (m !== 0 && (m < 1e-3 || m >= 1e5)) return v.toExponential(2);
    return String(Number(v.toPrecision(4)));
  }
  /* Eksen işaretleri: 1-2-5 adımı. Aralık sıfırsa çağıran önce genişletir. */
  function niceTicks(lo, hi, count = 5) {
    const span = hi - lo;
    if (!(span > 0) || !Number.isFinite(span)) return [lo];
    const raw = span / count, mag = 10 ** Math.floor(Math.log10(raw));
    const stepSize = [1, 2, 5, 10].map((f) => f * mag).find((s) => s >= raw);
    const out = [];
    for (let t = Math.ceil(lo / stepSize) * stepSize; t <= hi + stepSize * 1e-9; t += stepSize) out.push(Math.abs(t) < stepSize * 1e-9 ? 0 : Number(t.toPrecision(12)));
    return out;
  }
  /* Uç koşulu etiketi: verilmeyen uçta ne olduğu probleme göre DOĞRU yazılır (2. mertebe sınır değer: doğal α u′ = 0;
     başlangıç değer problemi ve 1. mertebe: serbest uç). */
  function conditionLabel(bc, problem) {
    const v = fmt(bc.value);
    switch (bc.type) {
      case 'dirichlet': return `u = ${v}`;
      case 'neumann': return `u′ = ${v}`;
      case 'flux': return `α u′ = ${v}`;
      case 'both': return `u = ${v}, u′ = ${fmt(bc.slope)}`;
      default: return problem.equationOrderX === 2 && !problem.initialValueProblem ? 'natural: α u′ = 0' : 'free end';
    }
  }

  /* ---------- Cebir aşamalarının çizimi ----------
     Matrisler ısı haritası: mavi pozitif, amber negatif, opaklık |değer| ile; büyük matriste yalnız sıfır olmayan hücreler çizilir
     (bantlı yapı: ≤ N(2p+1) hücre). Hücre yeterince büyükse değer yazılır. Değerler çekirdeğin sayılarıdır, çizici hesap yapmaz. */
  const ALGEBRA = new Set(['ansatz', 'residual', 'weak', 'elementK', 'assembly', 'bc', 'solve']);
  const cellText = (v) => {
    if (v === 0) return '0';
    const m = Math.abs(v);
    return m >= 1e-3 && m < 1e4 ? String(Number(v.toPrecision(3))) : v.toExponential(1);
  };
  function algebraFrame(data, fs, opts, W, Hreq) {
    const { stage, local: k, selNode, selEl } = fs;
    const D = data.discretization, A = data.algebra, B = data.basis, P = data.problem;
    const { elementCount: n, nodeCount: N, polynomialDegree: p, coordinates: X, connectivity: C, elementBoundaries: EB } = D;
    const [x0, x1] = P.domain, order = P.equationOrderX, S = A.termSamples;
    const angular = !!P.angularX;
    let H = Hreq || 300, body = '', note = '';
    const legend = [], r1 = (v) => v.toFixed(1);
    const text = (x, y, s, anchor = 'middle', cls = '') => `<text x="${r1(x)}" y="${r1(y)}" text-anchor="${anchor}"${cls ? ` class="${cls}"` : ''}>${s}</text>`;
    const enter = (str) => (str ? `<g class="enter">${str}</g>` : '');
    const sweep = (str) => (str ? `<g class="sweep">${str}</g>` : '');
    const pad = { l: 52, r: 16, t: 26, b: 32 };
    /* Fonksiyon eksenleri (rezidüel/zayıf form/çözme aşamalarındaki eğriler için) */
    const frameAxes = (lo, hi, yTitle, hgt = H) => {
      const pw = W - pad.l - pad.r, ph = hgt - pad.t - pad.b;
      const sx = (x) => pad.l + (x - x0) / (x1 - x0) * pw, sy = (y) => pad.t + (hi - y) / (hi - lo) * ph;
      let g = '';
      for (const t of niceTicks(lo, hi, 5)) g += `<line class="grid" x1="${pad.l}" x2="${W - pad.r}" y1="${r1(sy(t))}" y2="${r1(sy(t))}"/>` + text(pad.l - 8, sy(t) + 4, fmt(t), 'end');
      if (lo <= 0 && hi >= 0) g += `<line class="axis" x1="${pad.l}" x2="${W - pad.r}" y1="${r1(sy(0))}" y2="${r1(sy(0))}"/>`;
      for (const t of xTicks(x0, x1, Math.max(2, Math.floor(pw / 90)), angular)) g += text(sx(t.value), hgt - 8, t.label);
      g += text(pad.l - 16, hgt - 8, angular ? 'rad →' : 'x →', 'end', 'axis-title') + text(pad.l, 14, `↑ ${yTitle}`, 'start', 'axis-title');
      for (const x of EB) g += `<line class="vline" x1="${r1(sx(x))}" x2="${r1(sx(x))}" y1="${pad.t}" y2="${hgt - pad.b}"/>`;
      return { g, sx, sy, ph, pw };
    };
    const curve = (pts, sx, sy, style) => {
      let d = '', pen = 'M';
      for (const [x, y] of pts) { if (y === null || !Number.isFinite(y)) { pen = 'M'; continue; } d += `${pen}${r1(sx(x))},${r1(sy(y))}`; pen = 'L'; }
      return d ? `<path d="${d}" fill="none" stroke-linejoin="round" stroke-linecap="round" ${style}/>` : '';
    };
    const range = (values, include = [0]) => {
      let lo = Infinity, hi = -Infinity;
      for (const v of values.concat(include)) if (Number.isFinite(v)) { lo = Math.min(lo, v); hi = Math.max(hi, v); }
      if (!Number.isFinite(lo)) { lo = -1; hi = 1; }
      if (hi - lo < 1e-300) { lo -= 1; hi += 1; }
      const m = (hi - lo) * 0.08;
      return [lo - m, hi + m];
    };
    const localPiece = (kk, e, kind = 'N') => B.elementX[e - 1].map((x, q) => [x, kind === 'N' ? B.N[kk][q] : B.dN[kk][q] * 2 / D.h]);

    /* ----- matris çizimi ----- */
    const heat = (val, nr, nc, ox, oy, cell, o = {}) => {
      let max = o.max || 1, fills = '', fresh = '', frameLines = '', labels = '', hits = '';
      const small = nr * nc <= 400;
      for (let i = 0; i < nr; i++) for (let j = 0; j < nc; j++) {
        const v = val(i, j), x = ox + j * cell, y = oy + i * cell;
        const isNew = o.isNew && o.isNew(i, j);
        let piece = '';
        if (v !== null && v !== undefined && v !== 0) {
          const op = 0.16 + 0.74 * Math.min(1, Math.abs(v) / max);
          piece += `<rect x="${r1(x)}" y="${r1(y)}" width="${r1(cell)}" height="${r1(cell)}" style="fill:${v > 0 ? 'var(--c-ref)' : 'var(--c-err)'}" opacity="${op.toFixed(2)}"/>`;
        } else if (small) piece += `<rect x="${r1(x)}" y="${r1(y)}" width="${r1(cell)}" height="${r1(cell)}" style="fill:none;stroke:var(--grid)" stroke-width="0.8"/>`;
        if (cell >= 30 && v !== null && v !== undefined && (v !== 0 || o.zeros)) piece += text(x + cell / 2, y + cell / 2 + 3.5, cellText(v), 'middle', 'sim-cell');   // yapısal sıfırlar yazılmaz: boş K "0 0 0" gibi görünmesin
        if (o.hit) piece += `<rect class="hit" data-cell="${i},${j}" x="${r1(x)}" y="${r1(y)}" width="${r1(cell)}" height="${r1(cell)}" fill="transparent"/>`;
        if (isNew) fresh += piece; else fills += piece;
      }
      frameLines = `<rect x="${r1(ox)}" y="${r1(oy)}" width="${r1(nc * cell)}" height="${r1(nr * cell)}" fill="none" style="stroke:var(--ink-2)" stroke-width="1.2"/>`;
      if (o.hi) { const [hi, hj] = o.hi; frameLines += `<rect x="${r1(ox + hj * cell)}" y="${r1(oy + hi * cell)}" width="${r1(cell)}" height="${r1(cell)}" fill="none" style="stroke:var(--ink)" stroke-width="2"/>`; }
      if (o.strikeRows) for (const i of o.strikeRows) frameLines += `<line x1="${r1(ox)}" x2="${r1(ox + nc * cell)}" y1="${r1(oy + (i + 0.5) * cell)}" y2="${r1(oy + (i + 0.5) * cell)}" style="stroke:var(--c-mesh)" stroke-width="${Math.max(1, Math.min(3, cell / 3)).toFixed(1)}" stroke-dasharray="4 3"/>`;
      if (o.strikeCols) for (const j of o.strikeCols) frameLines += `<line y1="${r1(oy)}" y2="${r1(oy + nr * cell)}" x1="${r1(ox + (j + 0.5) * cell)}" x2="${r1(ox + (j + 0.5) * cell)}" style="stroke:var(--c-mesh)" stroke-width="${Math.max(1, Math.min(3, cell / 3)).toFixed(1)}" stroke-dasharray="4 3"/>`;
      const every = Math.max(1, Math.ceil(13 / cell));
      if (o.rowLabel) for (let i = 0; i < nr; i += every) labels += text(ox - 5, oy + (i + 0.5) * cell + 3.5, o.rowLabel(i), 'end', 'sim-node');
      if (o.colLabel && cell >= 12) for (let j = 0; j < nc; j += Math.max(every, Math.ceil(20 / cell))) labels += text(ox + (j + 0.5) * cell, oy - 5, o.colLabel(j), 'middle', 'sim-node');
      return fills + enter(fresh) + frameLines + labels;
    };
    /* K a = F üç blok: matris, bilinmeyen sütunu, sağ taraf sütunu. cell sığana kadar küçülür. */
    const system = (nr, nc, Kv, aV, Fv, o = {}) => {
      const maxCell = o.maxCell || 46, availH = H - pad.t - 30;
      const fits = (cell) => { const vw = Math.max(12, Math.min(cell, 44)); return nc * cell + 2 * 28 + 2 * vw <= W - 40 && Math.max(nr, nc) * cell <= availH; };
      let cell = maxCell;
      while (cell > 0.8 && !fits(cell)) cell -= cell > 6 ? 0.5 : 0.1;
      const vw = Math.max(12, Math.min(cell, 44)), total = nc * cell + 2 * 28 + 2 * vw, ox = Math.max(30, (W - total) / 2), oy = pad.t + 14;
      const aox = ox + nc * cell + 28, fox = aox + vw + 28;
      let g = heat(Kv, nr, nc, ox, oy, cell, { max: o.maxK, rowLabel: o.rowLabel, colLabel: o.colLabel, isNew: o.isNewK, strikeRows: o.strikeRows, strikeCols: o.strikeCols, hi: o.hiK, hit: false });
      g += text(ox + nc * cell / 2, oy - 18, o.kTitle || 'K', 'middle', 'axis-title');
      g += text(aox + vw / 2, oy - 18, 'a', 'middle', 'axis-title') + text(fox + vw / 2, oy - 18, o.fTitle || 'F', 'middle', 'axis-title');
      const mid = oy + Math.min(nr, nc) * cell / 2 + 5;
      g += text(aox - 14, mid, '×', 'middle', 'sim-op') + text(fox - 14, mid, '=', 'middle', 'sim-op');
      // a sütunu
      let acol = '';
      for (let i = 0; i < (o.aRows || nc); i++) {
        const av = aV(i), y = oy + i * cell;
        acol += `<rect x="${r1(aox)}" y="${r1(y)}" width="${r1(vw)}" height="${r1(cell)}" style="fill:${av && av.fixed ? 'var(--c-mesh)' : 'none'};stroke:var(--grid)" ${av && av.fixed ? 'opacity="0.35"' : ''} stroke-width="0.8"/>`;
        if (cell >= 14 && av) acol += text(aox + vw / 2, y + cell / 2 + 3.5, av.label, 'middle', 'sim-cell');
      }
      g += acol + `<rect x="${r1(aox)}" y="${r1(oy)}" width="${r1(vw)}" height="${r1((o.aRows || nc) * cell)}" fill="none" style="stroke:var(--ink-2)" stroke-width="1.2"/>`;
      g += heat((i) => Fv(i), nr, 1, fox, oy, cell, { max: o.maxF, isNew: o.isNewF ? (i) => o.isNewF(i) : null }).replace(/<rect x="[^"]*" y="[^"]*" width="[^"]*" height="[^"]*" fill="none" style="stroke:var\(--grid\)" stroke-width="0.8"\/>/g, '');
      g = g.replace('</svg>', '');
      // sütun genişliği: heat 1 sütunu cell genişlikte çizer; vw ile aynı olsun diye kalın çerçeveyi vw ile çiz
      g += `<rect x="${r1(fox)}" y="${r1(oy)}" width="${r1(cell)}" height="${r1(nr * cell)}" fill="none" style="stroke:var(--ink-2)" stroke-width="1.2"/>`;
      return { g, cell, ox, oy, rightX: fox + cell, bottom: oy + Math.max(nr, nc) * cell };
    };
    const aV = (j) => {
      const fixed = A.fixed.find(([idx]) => idx === j);
      return fixed ? { fixed: true, label: cellText(fixed[1]) } : { fixed: false, label: `a${String(j + 1).split('').map((c) => '₀₁₂₃₄₅₆₇₈₉'[+c]).join('')}` };
    };
    const maxOf = (rows) => { let m = 0; for (const r of rows) for (const v of Array.isArray(r) ? r : [r]) m = Math.max(m, Math.abs(v)); return m || 1; };
    const sub = (i) => String(i).split('').map((c) => '₀₁₂₃₄₅₆₇₈₉'[+c]).join('');

    if (stage.id === 'ansatz') {
      const yLine = pad.t + (H - pad.t - pad.b) * 0.62;
      body += `<line x1="${r1(pad.l)}" x2="${r1(W - pad.r)}" y1="${yLine}" y2="${yLine}" style="stroke:var(--ink)" stroke-width="2.4"/>`;
      const sx = (x) => pad.l + (x - x0) / (x1 - x0) * (W - pad.l - pad.r);
      if (k >= 1 && N <= 65) {
        const sc = frameAxes(-0.1, 1.2, 'N_j(x)');
        body = sc.g.replace(/<line class="axis"[^>]*\/>/, '') + body;
        for (let e = 1; e <= n; e++) for (let kk = 0; kk <= p; kk++) body += curve(localPiece(kk, e), sc.sx, (y) => yLine - y * (yLine - pad.t - 10) / 1.1, 'style="stroke:var(--c-basis-soft)" stroke-width="1.2"');
      }
      let chips = '', last = -1e9;
      for (let j = 1; j <= N; j++) {
        const px = sx(X[j - 1]);
        chips += `<circle cx="${r1(px)}" cy="${yLine}" r="3.6" style="fill:${(j - 1) % p ? 'var(--surface)' : 'var(--ink)'};stroke:var(--ink)" stroke-width="1.5"/>`;
        if (px - last >= 26 || j === N || j === 1) { chips += text(px, yLine + 26, `a${sub(j)}`, 'middle', 'sim-fn'); last = px; }
      }
      body += sweep(chips);
      note = k === 0 ? `${N} unknown coefficients a_j, one per node; u_h = Σ a_j N_j is the trial solution`
        : 'The N_j span the trial space; the residual R = 𝓛(u_h) − f is not zero, so it is made orthogonal to chosen weight functions';
    } else if (stage.id === 'weak' || stage.id === 'residual') {
      const hi = stage.id === 'weak' && k === 1 && order === 2;
      const sup = B.support[selNode - 1];
      const ders = [];
      if (hi) for (const { element, local: kk } of sup) ders.push(...localPiece(kk - 1, element, 'dN').map((q) => q[1]));
      const [lo, up] = hi ? range(ders) : range(B.N.flat(), [0, 1]);
      const sc = frameAxes(lo, up, hi ? `N_${selNode}′(x)` : 'weights w_i = N_i(x)');
      body += sc.g;
      if (!hi) {
        let faint = '';
        for (let e = 1; e <= n; e++) { if (N > 65 && !sup.some((q) => Math.abs(q.element - e) <= 2)) continue; for (let kk = 0; kk <= p; kk++) faint += curve(localPiece(kk, e), sc.sx, sc.sy, 'style="stroke:var(--c-basis-soft)" stroke-width="1"'); }
        let bold = '';
        for (const { element, local: kk } of sup) bold += curve(localPiece(kk - 1, element), sc.sx, sc.sy, 'style="stroke:var(--c-basis)" stroke-width="2.8" class="glow"');
        body += faint + sweep(bold);
        legend.push({ label: `w_${selNode} = N_${selNode}`, color: 'var(--c-basis)', width: 2.8 }, { label: 'other weights', color: 'var(--c-basis-soft)', width: 1.4 });
        note = `Galerkin: the weight function w_i is the same shape function N_i as the trial space (node ${selNode} selected; click a node on the basis card to change)`;
      } else {
        let d = '';
        for (const { element, local: kk } of sup) d += curve(localPiece(kk - 1, element, 'dN'), sc.sx, sc.sy, 'style="stroke:var(--c-deriv)" stroke-width="2.8" class="glow"');
        body += sweep(d);
        legend.push({ label: `N_${selNode}′ (derivative of the weight)`, color: 'var(--c-deriv)', width: 2.8 });
        note = 'After integration by parts the equation needs only first derivatives N_i′ and u_h′; the second derivative of u_h is never required, and N_j′ jumps at element boundaries';
      }
    } else if (stage.id === 'elementK') {
      const Ke = A.elementK[selEl - 1], Fe = A.elementF[selEl - 1], { kPhase, cell: cc, fIndex } = fs, e0 = selEl - 1;
      const filledK = (i, j) => kPhase === 'F' || kPhase === 'final' || (kPhase === 'K' && (i < cc.i || (i === cc.i && j <= cc.j)));
      const filledF = (i) => kPhase === 'final' || (kPhase === 'F' && i <= fIndex);
      const wide = W >= 520;
      if (!wide) H = Hreq ? Hreq + 220 : 520;
      const cell = Math.min(46, wide ? (H - pad.t - 60) / (p + 1) : (W - 40 - 28 - 46) / (p + 1)), ox = 54, oy = pad.t + 22;
      body += text(ox + (p + 1) * cell / 2, oy - 12, `K⁽${selEl}⁾`, 'middle', 'axis-title') + text(ox + (p + 1) * cell + 28 + cell / 2, oy - 12, `F⁽${selEl}⁾`, 'middle', 'axis-title');
      body += heat((i, j) => (filledK(i, j) ? Ke[i][j] : null), p + 1, p + 1, ox, oy, cell,
        { max: maxOf(Ke), hit: true, zeros: true, hi: kPhase === 'K' ? [cc.i, cc.j] : null, rowLabel: (i) => String(i + 1), colLabel: (j) => String(j + 1), isNew: (i, j) => kPhase === 'K' && i === cc.i && j === cc.j });
      body += heat((i) => (filledF(i) ? Fe[i] : null), p + 1, 1, ox + (p + 1) * cell + 28, oy, cell,
        { max: maxOf(Fe), zeros: true, hi: kPhase === 'F' ? [fIndex, 0] : null, isNew: (i) => kPhase === 'F' && i === fIndex });
      // sağ/alt grafik: seçili girdinin integrandı ve altındaki alan; giriş ve son adımda yerel şekil fonksiyonları
      const px0 = wide ? ox + (p + 1) * cell + 28 + cell + 40 : pad.l, py0 = wide ? pad.t : oy + (p + 1) * cell + 40;
      const pw = W - px0 - pad.r, ph = (wide ? H - pad.t - pad.b : H - py0 - pad.b);
      const xs = B.elementX[e0];
      const ix = (x) => px0 + (x - xs[0]) / (xs[xs.length - 1] - xs[0]) * pw;
      let series = [], areaVals = null, title = '', valueText = '';
      const col = ['--c-basis', '--c-deriv', '--c-err', '--c-ref', '--c-fem'];
      if (kPhase === 'K') {
        const q = (t) => xs.map((_, s2) => termAt(A, p, e0, cc.i, cc.j, t, s2)), tot = xs.map((_, s2) => q(0)[s2] + q(1)[s2] + q(2)[s2]);
        areaVals = tot; title = `integrand of K⁽${selEl}⁾${sub(cc.i + 1)}${sub(cc.j + 1)}`;
        [['--c-deriv', "N_i′ α N_j′", q(0)], ['--c-basis', '− N_i β N_j′', q(1)], ['--c-err', '− N_i γ N_j', q(2)]].forEach(([c, lab, v]) => {
          if (v.some((z) => Math.abs(z) > 0)) { series.push({ vals: v, color: c, width: 1.3, dash: '4 3' }); legend.push({ label: lab, color: `var(${c})`, width: 1.6, dash: '4 3' }); }
        });
        series.push({ vals: tot, color: '--c-deriv', width: 2.8, glow: true }); legend.push({ label: 'integrand (sum)', color: 'var(--c-deriv)', width: 2.8 });
        valueText = `K⁽${selEl}⁾${sub(cc.i + 1)}${sub(cc.j + 1)} = ${cellText(Ke[cc.i][cc.j])}`;
      } else if (kPhase === 'F') {
        const v = xs.map((_, s2) => A.fterms[(e0 * (p + 1) + fIndex) * S + s2]);
        areaVals = v; title = `integrand of F⁽${selEl}⁾${sub(fIndex + 1)}`;
        series.push({ vals: v, color: '--c-input', width: 2.8, glow: true }); legend.push({ label: order === 2 ? '− N_i f' : 'N_i f', color: 'var(--c-input)', width: 2.8 });
        valueText = `F⁽${selEl}⁾${sub(fIndex + 1)} = ${cellText(Fe[fIndex])}`;
      } else {
        title = `shape functions of element ${selEl}`;
        for (let kk = 0; kk <= p; kk++) { series.push({ vals: B.N[kk], color: col[kk % col.length], width: 2.2 }); legend.push({ label: `N${sub(kk + 1)}`, color: `var(${col[kk % col.length]})`, width: 2.2 }); }
      }
      const all = [0, ...series.flatMap((z) => z.vals), ...(areaVals || [])];
      const [lo, hi2] = range(all), iy = (y) => py0 + (hi2 - y) / (hi2 - lo) * ph;
      let plot = `<rect x="${r1(px0)}" y="${r1(py0)}" width="${r1(pw)}" height="${r1(ph)}" fill="none" style="stroke:var(--grid)"/>`;
      for (const t of niceTicks(lo, hi2, 3)) plot += `<line class="grid" x1="${r1(px0)}" x2="${r1(px0 + pw)}" y1="${r1(iy(t))}" y2="${r1(iy(t))}"/>` + text(px0 - 5, iy(t) + 4, fmt(t), 'end', 'sim-node');
      if (lo <= 0 && hi2 >= 0) plot += `<line class="axis" x1="${r1(px0)}" x2="${r1(px0 + pw)}" y1="${r1(iy(0))}" y2="${r1(iy(0))}"/>`;
      plot += text(px0 + pw, py0 + ph + 14, fmt(xs[xs.length - 1]), 'end', 'sim-node') + text(px0, py0 + ph + 14, fmt(xs[0]), 'start', 'sim-node') + text(px0, py0 - 7, title, 'start', 'axis-title');
      if (areaVals) {
        let area = `M${r1(ix(xs[0]))},${r1(iy(0))}`;
        xs.forEach((x, q2) => { area += `L${r1(ix(x))},${r1(iy(areaVals[q2]))}`; });
        plot += `<path d="${area}L${r1(ix(xs[xs.length - 1]))},${r1(iy(0))}Z" style="fill:var(--c-deriv)" opacity="0.16"/>`;
      }
      series.forEach((z, idx) => {
        const main = idx === series.length - 1 && areaVals;
        const c = curve(xs.map((x, q2) => [x, z.vals[q2]]), ix, iy, `style="stroke:var(${z.color});color:var(${z.color})" stroke-width="${z.width}"${z.dash ? ` stroke-dasharray="${z.dash}"` : ''}${z.glow ? ' class="glow"' : ''}`);
        plot += main || !areaVals ? sweep(c) : c;
      });
      body += plot;
      note = kPhase === 'intro' ? `Element ${selEl}: entry (i, j) of K⁽ᵉ⁾ pairs the weight N_i with the unknown N_j; click a cell to jump to it`
        : kPhase === 'final' ? `K⁽${selEl}⁾ and F⁽${selEl}⁾ complete`
          : `${valueText}: the shaded area under the integrand over element ${selEl}`;
    } else if (stage.id === 'assembly') {
      const upto = k, Mk = Array.from({ length: N }, () => new Array(N).fill(0)), Fk = new Array(N).fill(0);
      for (let e = 0; e < upto; e++) for (let i = 0; i <= p; i++) { Fk[e * p + i] += A.elementF[e][i]; for (let j = 0; j <= p; j++) Mk[e * p + i][e * p + j] += A.elementK[e][i][j]; }
      const inNew = (i, j) => upto > 0 && i >= (upto - 1) * p && i <= upto * p && j >= (upto - 1) * p && j <= upto * p;
      const sys = system(N, N, (i, j) => Mk[i][j], (j) => ({ fixed: false, label: `a${sub(j + 1)}` }), (i) => Fk[i],   // sınır koşulu henüz uygulanmadı: a_j hep simge
        { maxK: maxOf(A.K), maxF: maxOf(A.Fl), rowLabel: (i) => String(i + 1), colLabel: (j) => String(j + 1), isNewK: inNew, hiK: upto >= 2 ? [p, p] : null, isNewF: (i) => upto > 0 && i >= (upto - 1) * p && i <= upto * p, kTitle: 'K', fTitle: 'F_l' });
      H = Math.max(H, sys.bottom + 22);
      body += sys.g;
      note = k === 0 ? `The global K starts empty (every entry 0): N × N with N = n·p + 1 = ${N}`
        : `Element ${k} added: its ${p + 1} × ${p + 1} block lands at rows and columns ${(k - 1) * p + 1}–${k * p + 1}; the node shared with element ${Math.max(1, k - 1)} accumulates both contributions${k >= 2 ? ' (outlined entry: K₍ₚ₊₁,ₚ₊₁₎ = K⁽¹⁾ + K⁽²⁾)' : ''}${N > 41 ? ' · only nonzero entries are drawn' : ''}`;
    } else if (stage.id === 'bc') {
      const hasFlux = order === 2, kFlux = hasFlux ? 1 : -1, kEss = hasFlux ? 2 : 1, kRed = hasFlux ? 3 : 2;
      const useF = k >= kFlux && hasFlux ? A.F : A.Fl;
      const strikeRows = k >= kEss && k < kRed ? [...Array(N).keys()].filter((i) => !A.rows.includes(i)) : [];
      const strikeCols = k >= kEss && k < kRed ? A.fixed.map(([i]) => i) : [];
      if (k >= kRed) {
        const nr = A.rows.length, nc = A.free.length;
        const sys = system(nr, nc, (i, j) => A.reducedK[i][j], (j) => ({ fixed: false, label: aV(A.free[j]).label }), (i) => A.reducedF[i],
          { rowLabel: (i) => String(A.rows[i] + 1), colLabel: (j) => String(A.free[j] + 1), aRows: nc, kTitle: 'K_red', fTitle: 'F_red', maxK: maxOf(A.reducedK), maxF: maxOf(A.reducedF) });
        H = Math.max(H, sys.bottom + 22);
        body += enter(sys.g);
        note = `Reduced system: ${nr} equations in ${nc} unknowns; the known value${A.fixed.length === 1 ? '' : 's'} moved to the right-hand side`;
      } else {
        const sys = system(N, N, (i, j) => A.K[i][j], (j) => aV(j), (i) => useF[i],
          { rowLabel: (i) => String(i + 1), colLabel: (j) => String(j + 1), maxK: maxOf(A.K), maxF: maxOf(A.F), strikeRows, strikeCols, fTitle: k >= kFlux && hasFlux ? 'F' : 'F_l', isNewF: (i) => hasFlux && k === kFlux && A.Fb[i] !== 0 });
        H = Math.max(H, sys.bottom + 22);
        body += sys.g;
        if (k < kFlux || !hasFlux) note = k === 0 ? 'Assembled system before any boundary condition: K is singular until the essential conditions fix the rigid-body freedom' : '';
        if (hasFlux && k === kFlux) note = 'Natural (flux) conditions add F_b = [−Q₀, 0, …, 0, Q₁]ᵀ to the end entries of F; with no condition given, Q = 0';
        if (k === kEss) note = `Essential condition${A.fixed.length === 1 ? '' : 's'}: ${A.fixed.map(([i, v]) => `a${sub(i + 1)} = ${cellText(v)}`).join(', ') || 'none'} fix${A.fixed.length === 1 ? 'es' : ''} the coefficient${A.fixed.length === 1 ? '' : 's'}; struck rows and columns are removed`
          + (A.rows.length < N - A.fixed.length ? '' : '');
      }
    } else if (stage.id === 'solve') {
      if (k === 0) {
        const nr = A.rows.length, nc = A.free.length;
        const sys = system(nr, nc, (i, j) => A.reducedK[i][j], (j) => ({ fixed: false, label: aV(A.free[j]).label }), (i) => A.reducedF[i],
          { rowLabel: (i) => String(A.rows[i] + 1), colLabel: (j) => String(A.free[j] + 1), aRows: nc, kTitle: 'K_red', fTitle: 'F_red', maxK: maxOf(A.reducedK), maxF: maxOf(A.reducedF) });
        H = Math.max(H, sys.bottom + 22);
        body += sys.g;
        note = `${nr} × ${nc} system` + (A.condition ? ` · condition estimate κ₁ ≈ ${A.condition.toExponential(1)} (a lower-bound estimate)` : '');
      } else {
        const [lo, up] = range(A.a);
        const sc = frameAxes(lo, up, 'a_j');
        body += sc.g;
        let stems = '';
        A.a.forEach((v, j) => {
          const px = sc.sx(X[j]), fixedNode = A.fixed.some(([i]) => i === j);
          stems += `<line x1="${r1(px)}" x2="${r1(px)}" y1="${r1(sc.sy(0))}" y2="${r1(sc.sy(v))}" style="stroke:var(--c-fem)" stroke-width="${N > 65 ? 1.2 : 2.2}"/>`
            + `<circle cx="${r1(px)}" cy="${r1(sc.sy(v))}" r="${N > 65 ? 1.8 : 3.8}" style="fill:${fixedNode ? 'var(--surface)' : 'var(--c-fem)'};stroke:var(--c-fem)" stroke-width="1.6"/>`;
        });
        body += sweep(stems);
        legend.push({ label: 'a_j solved', color: 'var(--c-fem)', width: 2.2 }, { label: 'a_j given (essential)', color: 'var(--c-fem)', width: 1.2, dash: '2 2' });
        note = `Solved coefficients: ${N <= 9 ? A.a.map((v, j) => `a${sub(j + 1)} = ${cellText(v)}`).join(', ') : `${N} values`}; each a_j is the value u_h(x_j) at its node`;
      }
    }
    const svg = `<svg viewBox="0 0 ${W} ${H}" width="${W}" height="${H}" role="img" data-xranges="${pad.l},${W - pad.r}" aria-label="${stage.title}: ${note.replace(/"/g, '')}">${body}</svg>`;
    const picked = sheetPick(opts.sheet, stage, fs, { order, p, n });
    return { svg, note, legend, probe: null, stageIndex: fs.index, local: k, formulas: picked ? picked.tex : stageFormulas(data, fs, opts), formulaSource: picked ? picked.source : null, scale: { ranges: [[pad.l, W - pad.r]], x0, x1 } };
  }

  /* ---------- 4. ÇİZİM ----------
     Neden kendi küçük SVG çizicisi: ağ numaraları, eleman seçimi, referans eleman eşlemesi ve destek gölgesi gibi
     öğeler genel grafik fonksiyonunda (drawChart) yok; onu değiştirmek diğer grafikleri riske atardı. Renkler Stage B
     tokenlarıdır (--c-*). Dönen scale, tıklanan pikseli x'e çevirmek içindir. */
  function renderFrame(data, stages, step, opts = {}) {
    const W = Math.max(260, Math.floor(opts.width || 640)), H = opts.height || 300;
    const fs = frameState(data, stages, step, opts), { index, stage, local: k, selNode, selEl } = fs;
    if (ALGEBRA.has(stage.id) && data.algebra) return algebraFrame(data, fs, opts, W, opts.height || 300);
    const D = data.discretization, S = data.solution, B = data.basis;
    const { elementCount: n, nodeCount, polynomialDegree: p, coordinates: X, connectivity: C, elementBoundaries: EB } = D;
    const [x0, x1] = data.problem.domain;
    const showMesh = opts.showMesh !== false, showLabels = opts.showLabels !== false;
    const pad = { l: 52, r: 16, t: 26, b: 32 }, plotW = W - pad.l - pad.r, plotH = H - pad.t - pad.b;
    const sx = (x) => pad.l + (x - x0) / (x1 - x0) * plotW;
    const r1 = (v) => v.toFixed(1);
    let back = '', body = '', front = '', note = '';
    const legend = [];
    /* Akış: bu adımda YENİ gelen öğeler soldan sağa çizilir (sweep) ya da belirir (enter); önceden görünenler durağandır.
       Yalnız görünme biçimidir: çizilen değerler her zaman hesaplanmış değerlerdir, ara durum uydurulmaz. */
    let xRanges = [[pad.l, W - pad.r]];   // tıklanan pikseli x'e çevirmek için: [sol, sağ] piksel aralıkları (girdi panellerinde birden çok)
    const sweep = (str) => (str ? `<g class="sweep">${str}</g>` : '');
    const enter = (str) => (str ? `<g class="enter">${str}</g>` : '');
    const path = (pts, sy, style) => {
      let d = '', pen = 'M';
      for (const [x, y] of pts) {
        if (y === null || !Number.isFinite(y)) { pen = 'M'; continue; }   // eksik değer: çizgi kırılır, 0 sanılmaz
        d += `${pen}${r1(sx(x))},${r1(sy(y))}`; pen = 'L';
      }
      /* Ana eğriler (kalınlık ≥ 2.4) ışıltılı: rengi currentColor olarak da verilir ki hale aynı renkte olsun. */
      const main = /stroke-width="(\d+(\.\d+)?)"/.exec(style), glow = main && Number(main[1]) >= 2.4 && !/stroke-dasharray="5 4"/.test(style);
      const col = /stroke:([^;"]+)/.exec(style);
      return d ? `<path d="${d}" fill="none"${glow ? ' class="glow"' : ''} stroke-linejoin="round" stroke-linecap="round" ${col ? style.replace(/style="stroke:([^;"]+)/, `style="stroke:$1;color:$1`) : style}/>` : '';
    };
    const text = (x, y, s, anchor = 'middle', cls = '') => `<text x="${r1(x)}" y="${r1(y)}" text-anchor="${anchor}"${cls ? ` class="${cls}"` : ''}>${s}</text>`;
    const nodeDot = (cx, cy, j, hollow, color, r = 3.6) =>
      `<circle cx="${r1(cx)}" cy="${r1(cy)}" r="${r}" style="fill:${hollow ? 'var(--surface)' : color};stroke:${color}" stroke-width="1.6"/>`
      + `<circle class="hit" data-node="${j}" cx="${r1(cx)}" cy="${r1(cy)}" r="9" fill="transparent"><title>node ${j}: x = ${fmt(X[j - 1])}</title></circle>`;
    /* Etiket seyreltme: öncelik sırasıyla (seçili, ilk, son, sonra soldan sağa) bir etiket, kabul edilmiş her etiketten en az
       minPx uzaksa yazılır; minPx en uzun etiketin rakam sayısından gelir. Böylece etiketler hiçbir zaman üst üste binmez. */
    const labelSet = (items, minPx) => {
      const taken = [], out = new Set();
      for (const { id, px } of items) if (!out.has(id) && taken.every((q) => Math.abs(q - px) >= minPx)) { out.add(id); taken.push(px); }
      return out;
    };
    const digitPx = (k) => 7 * String(k).length + 8;   // 10.5 px eş aralıklı yazıda bir rakam ≈ 6.3 px
    const angular = !!data.problem.angularX;
    const xTicksSvg = (y) => xTicks(x0, x1, Math.max(2, Math.floor(plotW / 90)), angular).map((t) => text(sx(t.value), y, t.label)).join('') + text(pad.l - 16, y, angular ? 'rad →' : 'x →', 'end', 'axis-title');
    const meshLines = (top, bottom) => (showMesh ? EB.map((x) => `<line class="vline" x1="${r1(sx(x))}" x2="${r1(sx(x))}" y1="${top}" y2="${bottom}"/>`).join('') : '');
    const hitElements = (top, bottom) => EB.slice(0, -1).map((xa, e) =>
      `<rect class="hit" data-element="${e + 1}" x="${r1(sx(xa))}" y="${top}" width="${r1(sx(EB[e + 1]) - sx(xa))}" height="${bottom - top}" fill="transparent"><title>element ${e + 1}</title></rect>`).join('');
    const elementShade = (e, top, bottom, color, opacity) =>
      `<rect x="${r1(sx(EB[e - 1]))}" y="${top}" width="${r1(sx(EB[e]) - sx(EB[e - 1]))}" height="${bottom - top}" style="fill:${color}" opacity="${opacity}"/>`;

    /* Dikey ölçek: valueRange (3B ile ortak) → piksel. */
    const yScale = () => {
      const { lo, hi } = valueRange(data, stage.id);
      return { lo, hi, sy: (y) => pad.t + (hi - y) / (hi - lo) * plotH };
    };
    const yAxis = (sc, title) => {
      let s = '';
      for (const t of niceTicks(sc.lo, sc.hi, 5)) s += `<line class="grid" x1="${pad.l}" x2="${W - pad.r}" y1="${r1(sc.sy(t))}" y2="${r1(sc.sy(t))}"/>` + text(pad.l - 8, sc.sy(t) + 4, fmt(t), 'end');
      if (sc.lo <= 0 && sc.hi >= 0) s += `<line class="axis" x1="${pad.l}" x2="${W - pad.r}" y1="${r1(sc.sy(0))}" y2="${r1(sc.sy(0))}"/>`;
      return s + text(pad.l, 14, `↑ ${title}`, 'start', 'axis-title');
    };
    const localN = (kk, e) => B.elementX[e - 1].map((x, s) => [x, B.N[kk][s]]);   // N_kk (0 tabanlı) e elemanında
    /* Girdi fonksiyonları: simülasyon yazılımlarındaki "problem data" görünümü gibi küçük çoklu grafik. Her fonksiyonun KENDİ
       dikey ölçeği vardır (birimleri ve büyüklükleri farklı; aynı eksen yanıltırdı). Sabitler yatay çizgi + değer. */
    const inputPanels = () => {
      const list = data.problem.inputFunctions || [], m = list.length;
      if (!m) return '';
      const cols = plotW >= 520 ? m : Math.min(2, m), rows = Math.ceil(m / cols), gap = 14;
      xRanges = [];
      const left = 18, top = 8, width = W - left - pad.r, height = H - top - 10;
      const pw = (width - gap * (cols - 1)) / cols, ph = (height - gap * (rows - 1)) / rows;
      let out = '';
      list.forEach((fn, i) => {
        const cx = left + (i % cols) * (pw + gap), cy = top + Math.floor(i / cols) * (ph + gap);
        const ix0 = cx + 40, ix1 = cx + pw - 6, iy0 = cy + 22, iy1 = cy + ph - 18;
        let lo = Infinity, hi = -Infinity;
        for (const v of fn.values) if (v !== null) { lo = Math.min(lo, v); hi = Math.max(hi, v); }
        if (!Number.isFinite(lo)) { lo = -1; hi = 1; }
        if (hi - lo < 1e-12 * Math.max(1, Math.abs(hi))) { const c = Math.max(Math.abs(hi), 1) * 0.5; lo -= c; hi += c; }
        const mg = (hi - lo) * 0.1; lo -= mg; hi += mg;
        const px = (x) => ix0 + (x - x0) / (x1 - x0) * (ix1 - ix0), py = (v) => iy0 + (hi - v) / (hi - lo) * (iy1 - iy0);
        xRanges.push([ix0, ix1]);
        let panel = `<rect x="${r1(cx)}" y="${r1(cy)}" width="${r1(pw)}" height="${r1(ph)}" rx="6" style="fill:none;stroke:var(--grid)"/>`;
        for (const t of niceTicks(lo, hi, 3)) panel += `<line class="grid" x1="${r1(ix0)}" x2="${r1(ix1)}" y1="${r1(py(t))}" y2="${r1(py(t))}"/>` + text(ix0 - 6, py(t) + 4, fmt(t), 'end');
        if (lo <= 0 && hi >= 0) panel += `<line class="axis" x1="${r1(ix0)}" x2="${r1(ix1)}" y1="${r1(py(0))}" y2="${r1(py(0))}"/>`;
        /* x işaretleri: ana eksenle aynı sistem (radyanda π katları). Panel dar olduğu için en çok 3 işaret; uç etiketleri kenara hizalı. */
        for (const t of xTicks(x0, x1, Math.max(2, Math.floor((ix1 - ix0) / 70)), angular)) {
          const tx = px(t.value);
          panel += text(tx, iy1 + 14, t.label, tx - ix0 < 14 ? 'start' : ix1 - tx < 14 ? 'end' : 'middle');
        }
        panel += text(cx + 10, cy + 15, fn.constant !== null ? `${fn.name} = ${fmt(fn.constant)} (constant)` : fn.name, 'start', 'sim-fn');
        let d = '', pen = 'M';
        S.samples.x.forEach((x, s) => { const v = fn.values[s]; if (v === null) { pen = 'M'; return; } d += `${pen}${r1(px(x))},${r1(py(v))}`; pen = 'L'; });
        panel += d ? sweep(`<path d="${d}" fill="none" class="glow" stroke-linejoin="round" stroke-linecap="round" style="stroke:var(--c-input);color:var(--c-input)" stroke-width="2.6"/>`) : '';
        if (probeX !== null) panel += `<line x1="${r1(px(probeX))}" x2="${r1(px(probeX))}" y1="${r1(iy0)}" y2="${r1(iy1)}" style="stroke:var(--ink-2)" stroke-width="1"/>`;
        out += panel;
      });
      note = `Input functions as entered (${list.map((f) => f.name).join(', ')}), each on its own scale${angular ? '; x in radians' : ''}; the solver uses the same parsed expressions. Click to read values.`;
      legend.push({ label: 'input function', color: 'var(--c-input)', width: 2.6 });
      return out;
    };
    const probeX = opts.probe !== null && opts.probe !== undefined && stage.id === 'domain' && k === 0 ? probeSample(data, opts.probe)?.x ?? null : null;

    if (stage.id === 'domain' || stage.id === 'mesh' || stage.id === 'element') {
      const yLine = pad.t + plotH * (stage.id === 'element' ? 0.36 : 0.5);
      back += xTicksSvg(H - 8);
      if (stage.id === 'domain' && k === 0) {   // ilk kare: girdiğin fonksiyonlar (panel açılır açılmaz görünür)
        back = '';
        body += inputPanels();
      } else if (stage.id === 'domain') {
        let line = `<line x1="${r1(sx(x0))}" x2="${r1(sx(x1))}" y1="${yLine}" y2="${yLine}" style="stroke:var(--ink)" stroke-width="3" stroke-linecap="round"/>`;
        for (const [x, lab, anchor] of [[x0, `x₀ = ${fmt(x0)}`, 'start'], [x1, `x₁ = ${fmt(x1)}`, 'end']])   // uç etiketleri içeri doğru: kenarda kesilmesin
          line += `<line x1="${r1(sx(x))}" x2="${r1(sx(x))}" y1="${yLine - 10}" y2="${yLine + 10}" style="stroke:var(--ink)" stroke-width="2"/>` + text(sx(x), yLine + 26, lab, anchor);
        line += text((sx(x0) + sx(x1)) / 2, yLine - 16, 'Ω = (x₀, x₁)', 'middle', 'axis-title');
        body += k === 1 ? sweep(line) : line;
        if (k >= 2) {
          const [L, R] = data.problem.boundaryConditions;
          body += enter(text(sx(x0), yLine - 40, conditionLabel(L, data.problem), 'start', 'sim-bc') + text(sx(x1), yLine - 40, conditionLabel(R, data.problem), 'end', 'sim-bc'));
          note = 'Boundary conditions at x₀ and x₁.';
        } else note = 'The continuous interval.';
      } else {
        const shown = fs.shownElements, finalMesh = fs.finalMesh;
        body += `<line x1="${r1(sx(x0))}" x2="${r1(sx(x1))}" y1="${yLine}" y2="${yLine}" style="stroke:var(--c-mesh)" stroke-width="1.2" stroke-dasharray="3 4"/>`;
        const prev = stage.id === 'mesh' ? fs.prevShown : n;        // önceki adımda görünen eleman sayısı
        body += `<line x1="${r1(sx(x0))}" x2="${r1(sx(EB[prev]))}" y1="${yLine}" y2="${yLine}" style="stroke:var(--ink)" stroke-width="2.4" stroke-linecap="round"/>`;
        if (shown > prev) body += sweep(`<line x1="${r1(sx(EB[prev]))}" x2="${r1(sx(EB[shown]))}" y1="${yLine}" y2="${yLine}" style="stroke:var(--ink)" stroke-width="2.4" stroke-linecap="round"/>`);
        let ticks = '', newTicks = '';
        for (let e = 0; e <= shown; e++) {
          const t = `<line x1="${r1(sx(EB[e]))}" x2="${r1(sx(EB[e]))}" y1="${yLine - 9}" y2="${yLine + 9}" style="stroke:var(--c-mesh)" stroke-width="1.2"/>`;
          if (e > prev) newTicks += t; else ticks += t;
        }
        body += ticks + enter(newTicks);
        const elMid = (e) => (sx(EB[e - 1]) + sx(EB[e])) / 2;
        const elOrder = [...(stage.id === 'element' ? [selEl] : []), 1, shown, ...Array.from({ length: shown }, (_, i) => i + 1)];
        const elLabels = labelSet(elOrder.map((e) => ({ id: e, px: elMid(e) })), digitPx(n) + 7);
        if (showLabels) for (const e of elLabels) {
          const lab = text(elMid(e), yLine + 24, `e${e}`, 'middle', stage.id === 'element' && e === selEl ? 'sim-sel' : 'sim-el');
          body += e > prev ? enter(lab) : lab;
        }
        if (stage.id === 'element') body = (k === 0 ? enter : (v) => v)(elementShade(selEl, pad.t, H - pad.b, 'var(--c-fem)', 0.1)) + body;
        const nodeOrder = [...(stage.id === 'element' ? C[selEl - 1] : []), 1, nodeCount, ...Array.from({ length: nodeCount }, (_, i) => i + 1)];
        const nodeLabels = labelSet(nodeOrder.map((j) => ({ id: j, px: sx(X[j - 1]) })), digitPx(nodeCount));
        for (let j = 1; j <= nodeCount; j++) {
          const end = (j - 1) % p === 0, inShown = X[j - 1] <= EB[shown] + 1e-12 * (x1 - x0);
          if (!inShown || (!end && !finalMesh)) continue;
          const inSel = stage.id === 'element' && C[selEl - 1].includes(j);
          /* Yeni: bu adımda açılan elemanın düğümü; son ağ adımında iç düğümler ve bütün numaralar. */
          const isNew = stage.id === 'mesh' && X[j - 1] > EB[prev] + 1e-12 * (x1 - x0);
          const lastMesh = stage.id === 'mesh' && k === stage.steps - 1;
          const dot = nodeDot(sx(X[j - 1]), yLine, j, !end, inSel ? 'var(--c-fem)' : 'var(--ink)', inSel ? 4.6 : 3.4);
          const lab = finalMesh && showLabels && nodeLabels.has(j) ? text(sx(X[j - 1]), yLine - 14, String(j), 'middle', inSel ? 'sim-sel' : 'sim-node') : '';
          front += (isNew || (lastMesh && !end) ? enter(dot) : dot) + (isNew || lastMesh ? enter(lab) : lab);
        }
        if (stage.id === 'mesh') {
          note = finalMesh ? `${n} element${n === 1 ? '' : 's'}, ${nodeCount} nodes: global node j = (e − 1)·p + k${p > 1 ? ' (hollow: interior nodes)' : ''}`
            + (showLabels && nodeLabels.size < nodeCount ? ` · ${nodeLabels.size} of ${nodeCount} node numbers shown (no overlap)` : '')
            : `Elements 1–${shown} of ${n}, h = ${fmt(D.h)}`;
        } else {
          /* Referans eleman s ∈ [0, h] (kâğıdın gösterimi) → seçili eleman: x = x_e + s. Yerel k → global (e − 1)p + k. İç çözücü ξ ∈ [−1, 1] kullanır: s = h(1 + ξ)/2. */
          const yRef = pad.t + plotH * 0.86, mid = (sx(EB[selEl - 1]) + sx(EB[selEl])) / 2;
          const half = Math.min(Math.max((sx(EB[selEl]) - sx(EB[selEl - 1])) / 2, 90, 30 * (p + 1)), plotW / 2);   // p + 1 etiket sığsın
          const ra = Math.max(pad.l, Math.min(mid - half, W - pad.r - 2 * half)), rb = ra + 2 * half;
          const rx = (t) => ra + (t + 1) / 2 * (rb - ra);
          if (k >= 1) {
            let inset = `<line x1="${r1(ra)}" x2="${r1(rb)}" y1="${yRef}" y2="${yRef}" style="stroke:var(--ink-2)" stroke-width="2"/>` + text(ra - 6, yRef + 4, 's', 'end', 'axis-title');
            let links = '';
            B.referenceNodes.forEach((t, kk) => {
              const j = C[selEl - 1][kk], px = sx(X[j - 1]);
              links += `<line x1="${r1(rx(t))}" x2="${r1(px)}" y1="${yRef - 4}" y2="${yLine + 5}" style="stroke:var(--c-mesh)" stroke-width="1" stroke-dasharray="2 3"/>`;
              inset += `<circle cx="${r1(rx(t))}" cy="${yRef}" r="3.4" style="fill:${kk === 0 || kk === p ? 'var(--ink-2)' : 'var(--surface)'};stroke:var(--ink-2)" stroke-width="1.4"/>`;
              if (showLabels) inset += text(rx(t), yRef + 18, `${fmt(D.h * (1 + t) / 2)}`, 'middle') + text(rx(t), yRef - 10, `${kk + 1}→${j}`, 'middle', 'sim-node');
            });
            body += enter(links);
            front += sweep(inset);
            note = `Element ${selEl}: x = ${fmt(EB[selEl - 1])} + s, 0 ≤ s ≤ h = ${fmt(D.h)} (Jacobian 1); labels k→j: local node k → global node j = (e − 1)·p + k`;
          } else note = `Element ${selEl} of ${n}: [${fmt(EB[selEl - 1])}, ${fmt(EB[selEl])}], local nodes ${C[selEl - 1].join(', ')} (click an element to select it)`;
        }
        back += hitElements(pad.t, H - pad.b);
      }
    } else if (stage.id === 'local' || stage.id === 'global') {
      const sc = yScale();
      back += yAxis(sc, stage.id === 'local' ? `N_k(x) on element ${selEl}` : 'N_j(x)') + xTicksSvg(H - 8) + meshLines(pad.t, H - pad.b);
      if (stage.id === 'local') {
        back += elementShade(selEl, pad.t, H - pad.b, 'var(--c-basis)', 0.07);
        for (let kk = 0; kk <= k; kk++) {
          const pth = path(localN(kk, selEl), sc.sy, `style="stroke:${kk === k ? 'var(--c-basis)' : 'var(--c-basis-soft)'}" stroke-width="${kk === k ? 2.8 : 1.6}"`);
          body += kk === k ? sweep(pth) : pth;
        }
        C[selEl - 1].forEach((j, kk) => { front += nodeDot(sx(X[j - 1]), sc.sy(kk === k ? 1 : 0), j, kk !== 0 && kk !== p, 'var(--c-basis)'); });
        note = `N_${k + 1} on element ${selEl}: equals 1 at local node ${k + 1} (global ${C[selEl - 1][k]}) and 0 at the other local nodes`;
      } else {
        const sup = B.support[selNode - 1];
        const neighbors = new Set(fs.neighbors);
        const many = n * (p + 1) > 400;   // çok parçada soluk zemin yalnız seçilen çevresinde çizilir
        let faint = '';
        for (let e = 1; e <= n; e++) {
          if (many && !sup.some((s) => Math.abs(s.element - e) <= 2)) continue;
          for (let kk = 0; kk <= p; kk++) faint += path(localN(kk, e), sc.sy, 'style="stroke:var(--c-basis-soft)" stroke-width="1"');
        }
        body += k === 0 ? sweep(faint) : faint;
        if (k >= 1) {
          let shade = '', bold = '';
          for (const { element } of sup) shade += elementShade(element, pad.t, H - pad.b, 'var(--c-basis)', 0.1);
          for (const { element, local: kk } of sup) bold += path(localN(kk - 1, element), sc.sy, 'style="stroke:var(--c-basis)" stroke-width="2.8"');
          back += k === 1 ? enter(shade) : shade;
          front += k === 1 ? sweep(bold) : bold;
        }
        if (k >= 2) {
          let nb = '';
          for (const j of neighbors) for (const { element, local: kk } of B.support[j - 1])
            nb += path(localN(kk - 1, element), sc.sy, 'style="stroke:var(--c-basis)" stroke-width="1.4" stroke-dasharray="5 4"');
          body += sweep(nb);
        }
        for (let j = 1; j <= nodeCount; j++) if (!many || neighbors.has(j) || j === selNode)
          front += nodeDot(sx(X[j - 1]), sc.sy(j === selNode && k >= 1 ? 1 : 0), j, (j - 1) % p !== 0, j === selNode ? 'var(--c-basis)' : 'var(--c-mesh)', j === selNode ? 4.4 : 2.8);
        const supText = `support: element${sup.length === 1 ? '' : 's'} ${sup.map((s) => s.element).join(', ')} = [${fmt(EB[sup[0].element - 1])}, ${fmt(EB[sup[sup.length - 1].element])}]`;
        note = k === 0 ? `All local shape functions of the ${n} element${n === 1 ? '' : 's'} (click a node to select N_j)`
          : k === 1 ? `N_${selNode} at x = ${fmt(X[selNode - 1])}: ${sup.map((s) => `N_${s.local} of element ${s.element}`).join(' + ')} · ${supText}`
            : `Neighbouring basis functions (dashed) overlap N_${selNode}: nodes ${[...neighbors].sort((a, b) => a - b).join(', ')}`;
        legend.push({ label: `N_${selNode} (selected)`, color: 'var(--c-basis)', width: 2.8 }, { label: 'other N_j', color: 'var(--c-basis-soft)', width: 1.4 });
      }
    } else if (stage.id === 'approx') {
      const a = S.nodalCoefficients;
      const { terms, complete } = fs;
      const sc = yScale();
      back += yAxis(sc, 'u_h(x)') + xTicksSvg(H - 8) + meshLines(pad.t, H - pad.b);
      /* a_j N_j katkıları (yeşil kesikli) ve kısmi toplam Σ_{j ≤ m} a_j N_j (turuncu): kısmi toplam gerçek bir ara büyüklüktür,
         "u_h" diye ETİKETLENMEZ; son adımda m = N olunca u_h'nin kendisidir. */
      const drawFrom = Math.max(1, terms - (nodeCount > 65 ? 64 : nodeCount) + 1);
      let oldParts = '', newParts = '';
      for (let j = drawFrom; j <= terms; j++) for (const { element, local: kk } of B.support[j - 1]) {
        const newest = !complete && j >= fs.newestFrom;
        const pth = path(B.elementX[element - 1].map((x, s) => [x, a[j - 1] * B.N[kk - 1][s]]), sc.sy,
          `style="stroke:${newest || j === selNode ? 'var(--c-basis)' : 'var(--c-basis-soft)'}" stroke-width="${newest || j === selNode ? 1.8 : 1}" stroke-dasharray="5 4"`);
        if (newest) newParts += pth; else oldParts += pth;
      }
      body += oldParts + sweep(newParts);
      const partial = (m) => {
        let out = '';
        for (let e = 1; e <= n; e++) {
          const pts = B.elementX[e - 1].map((x, s) => [x, C[e - 1].reduce((acc, j, kk) => acc + (j <= m ? a[j - 1] * B.N[kk][s] : 0), 0)]);
          out += path(pts, sc.sy, `style="stroke:var(--c-fem)" stroke-width="${complete ? 3.2 : 2.4}"`);
        }
        return out;
      };
      if (fs.prevTerms) body += `<g opacity="0.28">${partial(fs.prevTerms)}</g>`;   // bir önceki kısmi toplam: değişimin izi
      front += complete ? enter(partial(terms)) : sweep(partial(terms));   // son adım: tam u_h kalınlaşarak belirir
      let dots = '', newDots = '';
      for (let j = 1; j <= terms; j++) if (nodeCount <= 129 || j === selNode) {
        const d = nodeDot(sx(X[j - 1]), sc.sy(a[j - 1]), j, (j - 1) % p !== 0, 'var(--c-fem)', j === selNode ? 4.6 : 3);
        if (!complete && j >= fs.newestFrom) newDots += d; else dots += d;
      }
      front += dots + enter(newDots);
      note = complete ? `u_h complete: all ${nodeCount} terms; at node ${selNode}, u_h(x_${selNode}) = a_${selNode} = ${fmt(a[selNode - 1])}`
        : `Partial sum of ${terms} of ${nodeCount} terms: Σ_{j ≤ ${terms}} a_j N_j(x) (not yet u_h)`;
      legend.push({ label: complete ? 'u_h = Σ a_j N_j' : 'partial sum', color: 'var(--c-fem)', width: 3 }, { label: 'a_j N_j', color: 'var(--c-basis)', width: 1.8, dash: '5 4' });
      if (fs.prevTerms) legend.push({ label: `previous partial sum (${fs.prevTerms} terms)`, color: 'var(--c-fem)', width: 2, faint: true });
    } else if (stage.id === 'reference') {
      const ref = S.samples.reference, xs = S.samples.x, analytical = data.metadata.referenceType === 'analytical';
      const sc = yScale();
      back += yAxis(sc, 'u(x)') + xTicksSvg(H - 8) + meshLines(pad.t, H - pad.b);
      const femPath = path(xs.map((x, i) => [x, S.samples.fem[i]]), sc.sy, `style="stroke:var(--c-fem)" stroke-width="3.2"${k === 0 ? ' opacity="0.35"' : ''}`);
      body += k === 1 ? enter(femPath) : femPath;
      const refPath = path(xs.map((x, i) => [x, ref[i]]), sc.sy, 'style="stroke:var(--c-ref)" stroke-width="2.8" stroke-dasharray="10 6"');
      front += k === 0 ? sweep(refPath) : refPath;
      if (k >= 1) {
        let dots = '';
        for (let j = 1; j <= nodeCount; j++) if (nodeCount <= 129) dots += nodeDot(sx(X[j - 1]), sc.sy(S.nodalCoefficients[j - 1]), j, (j - 1) % p !== 0, 'var(--c-fem)', 2.8);
        front += enter(dots);
      }
      note = (analytical ? `analytical u(x) (${data.metadata.referenceLevel || 'verified'})` : 'fine-mesh reference u_ref(x)') + (k >= 1 ? ' and u_h(x)' : '');
      legend.push({ label: analytical ? 'analytical u' : 'reference u_ref', color: 'var(--c-ref)', width: 2.8, dash: '10 6' }, { label: 'FEM u_h', color: 'var(--c-fem)', width: 3.2 });
    } else if (stage.id === 'error') {
      const err = S.samples.error, xs = S.samples.x, exact = data.metadata.errorKind === 'exact', sym = exact ? 'e' : 'd';
      const sc = yScale();
      back += yAxis(sc, exact ? 'e(x) = u − u_h' : 'd(x) = u_ref − u_h') + xTicksSvg(H - 8) + meshLines(pad.t, H - pad.b);
      const errPath = path(xs.map((x, i) => [x, err[i]]), sc.sy, 'style="stroke:var(--c-err)" stroke-width="2.4"');
      body += k === 0 ? sweep(errPath) : errPath;
      let peak = 0;
      for (const v of err) if (v !== null) peak = Math.max(peak, Math.abs(v));
      if (k >= 1 && peak > 0) {
        /* Büyük hata bölgesi: örneklenmiş |e|'nin eleman içi en büyüğü, genel en büyüğün yarısını aşan elemanlar gölgelenir. */
        const big = [];
        for (let e = 1; e <= n; e++) {
          let m = 0;
          xs.forEach((x, i) => { if (x >= EB[e - 1] - 1e-15 && x <= EB[e] + 1e-15 && err[i] !== null) m = Math.max(m, Math.abs(err[i])); });
          if (m >= 0.5 * peak) { big.push(e); back += enter(elementShade(e, pad.t, H - pad.b, 'var(--c-err)', 0.12)); }
        }
        const me = data.metadata.maxError;
        if (me) front += enter(`<line x1="${r1(sx(me.x))}" x2="${r1(sx(me.x))}" y1="${pad.t}" y2="${H - pad.b}" style="stroke:var(--c-err)" stroke-width="1.2" stroke-dasharray="2 3"/>`);
        note = `shaded: elements where max |${sym}| ≥ half of the largest value (elements ${big.join(', ')})`
          + (me ? ` · max |${sym}| ≈ ${fmt(me.value)} at x ≈ ${fmt(me.x)} (numerical search; lower estimate)` : '');
      } else note = peak > 0 ? `${sym}(x) at ${xs.length} sample points` : `${sym}(x) = 0 at every sample point`;
      legend.push({ label: exact ? 'error e(x)' : 'difference d(x)', color: 'var(--c-err)', width: 2.4 });
    }

    /* Sorgu noktası: tıklanan x'e en yakın ÖRNEK (gerçek hesaplanmış değerler); fonksiyon aşamalarında çizilir. */
    let probe = null;
    if (opts.probe !== null && opts.probe !== undefined && ['approx', 'reference', 'error'].includes(stage.id)) {
      probe = probeSample(data, opts.probe);
      if (probe) front += `<line x1="${r1(sx(probe.x))}" x2="${r1(sx(probe.x))}" y1="${pad.t}" y2="${H - pad.b}" style="stroke:var(--ink-2)" stroke-width="1"/>`;
    } else if (probeX !== null) probe = probeSample(data, opts.probe);
    const svg = `<svg viewBox="0 0 ${W} ${H}" width="${W}" height="${H}" role="img" data-xranges="${xRanges.map(([l, r]) => `${r1(l)},${r1(r)}`).join(';')}" aria-label="${stage.title}: ${note.replace(/"/g, '')}">`
      + `${back}${body}${front}</svg>`;
    const picked = sheetPick(opts.sheet, stage, fs, { order: data.problem.equationOrderX, p, n });
    return { svg, note, legend, probe, stageIndex: index, local: k, formulas: picked ? picked.tex : stageFormulas(data, fs, opts), formulaSource: picked ? picked.source : null, scale: { ranges: xRanges, x0, x1 } };
  }

  /* ---------- 5. SORGU ----------
     Tıklanan x'e en yakın örnek noktadaki HESAPLANMIŞ değerler (yeni hesap yok): x, u_h, referans, hata/fark. */
  function probeSample(data, x) {
    const xs = data.solution.samples.x;
    if (!Number.isFinite(x)) return null;
    let best = 0;
    for (let i = 1; i < xs.length; i++) if (Math.abs(xs[i] - x) < Math.abs(xs[best] - x)) best = i;
    const S = data.solution.samples;
    return { index: best, x: xs[best], fem: S.fem[best], reference: S.reference ? S.reference[best] : null, error: S.error ? S.error[best] : null,
      inputs: (data.problem.inputFunctions || []).map((fn) => ({ name: fn.name, value: fn.values[best] })) };
  }

  return { buildSimulationData, buildStages, stagesFor, locateStep, frameState, valueRange, renderFrame, stageFormulas, parseSheet, sheetPick, probeSample, niceTicks, xTicks, fmt, conditionLabel };
}

if (typeof module !== 'undefined' && module.exports) module.exports = { femSimulation, usesTrigonometry, angularTicks, parseSheet, sheetPick };
