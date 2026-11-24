// basis_core.js testleri (ana çözücünün Basis functions kartı):  node tests_basis.js
// Beklenen değerler çekirdekten DEĞİL, kitabın kendi formüllerinden (aşağıda harfiyen yazılı) ve elle/analitik hesaptan gelir.
// Mevcut tests.js'e dokunulmaz; bu dosya ona eklenen taban fonksiyonu işlevlerini sınar.
const fs = require('fs');
const { femExtension, coreFromHtml } = require('./basis_core.js');
const core = coreFromHtml(fs.readFileSync(__dirname + '/index.html', 'utf8'));
const X = femExtension(core);

let failures = 0;
function check(name, ok, info = '') { console.log(`${ok ? 'OK  ' : 'FAIL'} ${name} ${info}`); if (!ok) failures++; }
const close = (a, b, tol = 1e-12) => Math.abs(a - b) <= tol * Math.max(1, Math.abs(b));
const matClose = (A, B, tol = 1e-12) => A.length === B.length && A.every((row, i) => row.every((v, j) => close(v, B[i][j], tol)));

// 1. Kitabın formülleri, harfiyen (bağımsız kaynak). Yerel x ∈ [0, h]: (3.14)–(3.15), (3.19), (3.25). ξ: (3.66), (3.75)–(3.77).
const BOOK_LOCAL = {
  1: [(x, h) => 1 - x / h, (x, h) => x / h],
  2: [(x, h) => 1 - 3 * x / h + 2 * (x / h) ** 2, (x, h) => 4 * (x / h) * (1 - x / h), (x, h) => (x / h) * (2 * x / h - 1)],
  3: [(x, h) => (1 - 3 * x / h) * (1 - 3 * x / (2 * h)) * (1 - x / h), (x, h) => (9 * x / h) * (1 - 3 * x / (2 * h)) * (1 - x / h),
      (x, h) => -(9 * x / (2 * h)) * (1 - 3 * x / h) * (1 - x / h), (x, h) => (x / h) * (1 - 3 * x / h) * (1 - 3 * x / (2 * h))]
};
const BOOK_NATURAL = {
  1: [(t) => (1 - t) / 2, (t) => (1 + t) / 2],
  2: [(t) => t * (t - 1) / 2, (t) => 1 - t * t, (t) => t * (t + 1) / 2]
};
const xis = Array.from({ length: 41 }, (_, s) => -1 + s / 20);
for (const p of [1, 2, 3]) {
  const h = 0.37;   // keyfi eleman boyu: x = h(1+ξ)/2
  const ok = xis.every((t) => core.shapeFunctions(p, t).every((v, i) => close(v, BOOK_LOCAL[p][i](h * (1 + t) / 2, h), 1e-13)));
  check(`P${p}: çözücünün N_i(ξ) = kitabın yerel N_i^(e)(x), x = h(1+ξ)/2 (${['', '3.14–3.15', '3.19', '3.25'][p]})`, ok);
  if (BOOK_NATURAL[p]) check(`P${p}: doğal koordinatta kitapla aynı (${p === 1 ? '3.66' : '3.75–3.77'})`,
    xis.every((t) => core.shapeFunctions(p, t).every((v, i) => close(v, BOOK_NATURAL[p][i](t), 1e-14))));
  // türev: kitabın yerel ifadesinin merkezi farkı ile dN/dx = (2/h) dN/dξ  (3.67)
  const eps = 1e-6, okd = xis.slice(1, -1).every((t) => core.shapeDerivatives(p, t).every((d, i) => {
    const x = h * (1 + t) / 2, fd = (BOOK_LOCAL[p][i](x + eps, h) - BOOK_LOCAL[p][i](x - eps, h)) / (2 * eps);
    return Math.abs(d * 2 / h - fd) < 1e-6 * Math.max(1, Math.abs(fd));
  }));
  check(`P${p}: dN/dx = (2/h) dN/dξ, kitabın ifadesinin türeviyle aynı (3.67)`, okd);
  const table = X.nodalTable(p);
  check(`P${p}: N_i(ξ_j) = δ_ij (3.26)`, table.every((row, j) => row.every((v, i) => close(v, i === j ? 1 : 0, 1e-14))));
  check(`P${p}: yerel düğümler soldan sağa eşit aralıklı (Şekil 3.5, 3.25: 0, h/3, 2h/3, h)`,
    JSON.stringify(X.localNodes(p).map((v) => +v.toFixed(12))) === JSON.stringify(Array.from({ length: p + 1 }, (_, k) => +(-1 + 2 * k / p).toFixed(12))));
}
// Kübik ξ biçimi (kitapta yok; (3.25)'ten türetildi) ekranda gösterilen ifadenin sayısal karşılığı
const CUBIC_XI = [(t) => -(1 - t) * (1 - 9 * t * t) / 16, (t) => 9 * (1 - t * t) * (1 - 3 * t) / 16, (t) => 9 * (1 - t * t) * (1 + 3 * t) / 16, (t) => -(1 + t) * (1 - 9 * t * t) / 16];
const CUBIC_DXI = [(t) => (1 + 18 * t - 27 * t * t) / 16, (t) => 9 * (-3 - 2 * t + 9 * t * t) / 16, (t) => 9 * (3 - 2 * t - 9 * t * t) / 16, (t) => (-1 + 18 * t + 27 * t * t) / 16];
check('P3: ekrandaki ξ biçimi ve türevi çözücüyle aynı', xis.every((t) => core.shapeFunctions(3, t).every((v, i) => close(v, CUBIC_XI[i](t), 1e-14))
  && core.shapeDerivatives(3, t).every((v, i) => close(v, CUBIC_DXI[i](t), 1e-13))));
// P4 açık ifadeleri (kitapta yok; çözücünün Lagrange tanımı) ekranda gösterilen biçimle çözücüyle aynı
const P4_XI = [(t) => t * (t - 1) * (4 * t * t - 1) / 6, (t) => -4 / 3 * t * (t * t - 1) * (2 * t - 1), (t) => (t * t - 1) * (4 * t * t - 1),
  (t) => -4 / 3 * t * (t * t - 1) * (2 * t + 1), (t) => t * (t + 1) * (4 * t * t - 1) / 6];
const P4_DXI = [(t) => (16 * t ** 3 - 12 * t * t - 2 * t + 1) / 6, (t) => -4 / 3 * (8 * t ** 3 - 3 * t * t - 4 * t + 1), (t) => 16 * t ** 3 - 10 * t,
  (t) => -4 / 3 * (8 * t ** 3 + 3 * t * t - 4 * t - 1), (t) => (16 * t ** 3 + 12 * t * t - 2 * t - 1) / 6];
check('P4: ekranda gösterilen açık ifade ve türevi çözücünün Lagrange tabanıyla aynı', xis.every((t) => core.shapeFunctions(4, t).every((v, i) => close(v, P4_XI[i](t), 1e-13))
  && core.shapeDerivatives(4, t).every((v, i) => close(v, P4_DXI[i](t), 1e-12))) && X.SOLVER_P4.natural.length === 5);
check('kitap ifadeleri yalnız P1–P3 (P4 kitapta yok, uydurulmaz)', JSON.stringify(Object.keys(X.BOOK)) === '["1","2","3"]');
check('P5 reddedilir (çözücü de 1–4)', (() => { try { X.localNodes(5); return false; } catch (e) { return /1–4/.test(e.message); } })());
for (const n of [1, 3, 8]) {
  const pr = X.basisProperties(X.makeGrid(2.5, n, 4, -1.2));
  check(`P4 (çözücünün Lagrange tabanı), n=${n}, x₀ = −1.2: δ_IJ, Σ N_I = 1, Σ dN_I/dx = 0`, pr.delta < 1e-13 && pr.unity < 1e-13 && pr.slope < 1e-10,
    `${pr.delta.toExponential(1)} ${pr.unity.toExponential(1)} ${pr.slope.toExponential(1)}`);
}
{ // x₀ ≠ 0: düğüm konumları ve u_h kuruluşu aralığa doğru taşınır (ana sayfadaki [x₀, x₁] için)
  const g = X.makeGrid(3, 4, 2, -1), poly = (x) => 2 - x + 0.5 * x * x, T = g.x.map(poly);
  check('x₀ = −1: düğümler −1 … 2, u_h derece ≤ p polinomu tam verir', close(g.x[0], -1) && close(g.x[g.count - 1], 2)
    && Array.from({ length: 31 }, (_, k) => -1 + 3 * k / 30).every((x) => close(X.interpolateAt(g, T, x), poly(x), 1e-12)));
  const b = X.globalBasis(g, 3);
  check('x₀ = −1: ortak uç düğümün desteği [x₁, x₅] = [−1, 0.5]', close(b.supportX[0], -1) && close(b.supportX[1], 0.5));
}

// 2. Fiziksel eleman: x(ξ) = ½(1−ξ)x_i + ½(1+ξ)x_{i+1} (3.65), J = h/2 (3.68), dN/dx = (2/h) dN/dξ
for (const p of [1, 2, 3]) {
  const s = X.physicalSamples(p, 2, 3.5, 21);
  check(`P${p}: x(ξ) kitabın (3.65) dönüşümü ile aynı`, s.xi.every((t, k) => close(s.x[k], 0.5 * (1 - t) * 2 + 0.5 * (1 + t) * 3.5, 1e-14)));
  check(`P${p}: J = h/2 = 0.75 her ξ'de`, s.J.every((j) => close(j, 0.75, 1e-14)));
  check(`P${p}: dN/dx = dN/dξ / J`, s.dNdx.every((row, i) => row.every((v, k) => close(v, s.dN[i][k] / 0.75, 1e-14))));
}

// 3. Ağ, bağlantı, global taban (Şekil 3.3, 3.6)
const g2 = X.makeGrid(1, 3, 2);
check('P2, 3 eleman: e_i = {x_{2i−1}, x_{2i}, x_{2i+1}} (3.20)', JSON.stringify(g2.connectivity) === JSON.stringify([[1, 2, 3], [3, 4, 5], [5, 6, 7]]));
check('düğüm sayısı n·p + 1', [1, 2, 3].every((p) => X.makeGrid(1, 4, p).count === 4 * p + 1));
const gb = X.globalBasis(g2, 3);
check('ortak uç düğümün global fonksiyonu iki elemanda tanımlı, desteği [x_1, x_5]', JSON.stringify(gb.support) === '[1,2]' && close(gb.supportX[1], 2 / 3, 1e-14));
check('iç düğümün (x_2) desteği yalnız kendi elemanı', JSON.stringify(X.globalBasis(g2, 2).support) === '[1]');
check('global fonksiyon eleman sınırında sürekli (iki yandan aynı değer)', [1, 2, 3].every((p) => {
  const g = X.makeGrid(2, 4, p);
  return Array.from({ length: g.count }, (_, k) => k + 1).every((I) => [1, 2, 3].every((e) => close(X.globalValue(g, I, e * g.h, 'left'), X.globalValue(g, I, e * g.h, 'right'), 1e-13)));
}));
for (const p of [1, 2, 3]) for (const n of [1, 2, 5, 8]) {
  const pr = X.basisProperties(X.makeGrid(1.7, n, p));
  check(`P${p}, n=${n}: δ_IJ, Σ N_I = 1, Σ dN_I/dx = 0 (ölçülen en büyük sapma)`, pr.delta < 1e-13 && pr.unity < 1e-13 && pr.slope < 1e-11,
    `${pr.delta.toExponential(1)} ${pr.unity.toExponential(1)} ${pr.slope.toExponential(1)}`);
}
// 4. Yaklaşık çözüm: düğümde T_i; derece ≤ p polinomu tam temsil (kitap madde 4, s. 39)
for (const p of [1, 2, 3]) {
  const g = X.makeGrid(2, 3, p), poly = (x) => 1 - 2 * x + (p >= 2 ? 0.5 * x * x : 0) + (p >= 3 ? -0.25 * x ** 3 : 0);
  const T = g.x.map(poly), xs = Array.from({ length: 33 }, (_, k) => 2 * k / 32);
  check(`P${p}: T(x) = Σ N_i T_i derece ≤ ${p} polinomu tam verir`, xs.every((x) => close(X.interpolateAt(g, T, x), poly(x), 1e-12)));
  const segs = X.interpolate(g, T, 9);
  check(`P${p}: katkıların toplamı T(x)`, segs.every((s) => s.x.every((_, k) => close(s.parts.reduce((a, q) => a + q.values[k], 0), s.T[k], 1e-13))));
}
check('yanlış sayıda düğüm değeri reddedilir', (() => { try { X.interpolate(X.makeGrid(1, 2, 1), [1, 2]); return false; } catch (e) { return /3 nodal values/.test(e.message); } })());

// 4b. Ana sayfa (index.html) entegrasyonu: kartın kurduğu u_h = Σ a_j N_j, çözücünün kendi evaluate'iyle aynı olmalı
{ const html = fs.readFileSync(__dirname + '/index.html', 'utf8'), src = html.split('/*CORE-START*/')[1].split('/*CORE-END*/')[0];
  const full = new Function(src + '; return { makeMesh, buildProblem, METHODS, evaluate };')();
  const pb = full.buildProblem({ alpha: '1', beta: '', gamma: '1', f: 'x', x0: -1, x1: 2, left: { type: 'dirichlet', value: '0' }, right: { type: 'neumann', value: '1' } });
  let worst = 0;
  for (const [n, p] of [[1, 1], [3, 2], [5, 3], [7, 4], [64, 4]]) {
    const mesh = full.makeMesh(-1, 2, n, p), c = full.METHODS.galerkin.solve(pb, mesh).c, g = X.makeGrid(3, n, p, -1);
    for (let k = 0; k <= 300; k++) { const x = -1 + 3 * k / 300; worst = Math.max(worst, Math.abs(X.interpolateAt(g, c, x) - full.evaluate(mesh, c, x))); }
  }
  check('ana sayfa kartı: u_h = Σ a_j N_j çözücünün evaluate(mesh, c, x) ile aynı (P1–P4, x₀ = −1)', worst < 1e-14, worst.toExponential(1));
  check('ana sayfa: eklenti yoksa kart gizli, çözücü etkilenmez (BASIS koruması ve tek ekleme noktası)',
    /const BASIS = typeof femExtension === 'function' \?/.test(html) && /card\.hidden = !BASIS \|\| !view/.test(html) && html.includes('<script src="basis_core.js"></script>'));
}

console.log(`\n${failures ? failures + ' BAŞARISIZ' : 'Tümü geçti'}`);
process.exit(failures ? 1 : 0);
