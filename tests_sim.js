// fem_sim.js testleri (FEM Simulation paneli, 2B):  node tests_sim.js
// Sınanan riskler: (1) simülasyon verisi çözücünün verisinden ayrılır mı, (2) eksik veri uydurulur mu (0 yerine null),
// (3) durağan problemde fiziksel zaman iddiası, (4) her aşamanın her karesi geçerli SVG mi (NaN/undefined/Infinity yok),
// (5) adım ↔ aşama eşlemesi, (6) etiketler üst üste biner mi, (7) ağır durumda (n = 64, p = 4) süre, (8) sayfa entegrasyonu.
// Beklenen değerler çözücünün kendi fonksiyonlarından (CORE evaluate, mesh.nodeX) ve basis_core.js'ten gelir.
const fs = require('fs');
const { femExtension } = require('./basis_core.js');
const { femSimulation } = require('./fem_sim.js');

const html = fs.readFileSync(__dirname + '/index.html', 'utf8');
const src = html.split('/*CORE-START*/')[1].split('/*CORE-END*/')[0];
const F = new Function(src + '; return { makeMesh, buildProblem, METHODS, evaluate, sampleXs, shapeFunctions, shapeDerivatives, gaussLegendre, solveLinear, luSolve, elementSystem };')();
const BASIS = femExtension(F);
const SIM = femSimulation({ basis: BASIS, evaluate: F.evaluate });
const SIMA = femSimulation({ basis: BASIS, evaluate: F.evaluate, elementSystem: F.elementSystem });   // cebir aşamalarıyla (sayfadaki kurulum)
const { spawnSync } = require('child_process');

let failures = 0;
function check(name, ok, info = '') { console.log(`${ok ? 'OK  ' : 'FAIL'} ${name} ${info}`); if (!ok) failures++; }
const D = (type, value = '', slope = '') => ({ type, value, slope });

// Farklı denklem sınıfları: 2. mertebe Dirichlet, Neumann, akı, başlangıç değer, 1. mertebe; değişken katsayı, farklı aralık.
const PROBLEMS = [
  { name: "u'' + u = x, Dirichlet", inputs: { alpha: '1', beta: '', gamma: '1', f: 'x', x0: 0, x1: 1, left: D('dirichlet', '0'), right: D('dirichlet', '0') },
    exact: (x) => x - Math.sin(x) / Math.sin(1) },
  { name: "−u'' = sin(πx), Dirichlet + Neumann", inputs: { alpha: '-1', beta: '', gamma: '', f: 'sin(pi*x)', x0: 0, x1: 1, left: D('dirichlet', '0'), right: D('neumann', '1') } },
  { name: "((1+x²)u')' + exp(x) u = 1, akı", inputs: { alpha: '1+x^2', beta: '', gamma: 'exp(x)', f: '1', x0: -1, x1: 2, left: D('dirichlet', '1'), right: D('flux', '2') } },
  { name: "u' = cos(x), 1. mertebe", inputs: { alpha: '', beta: '1', gamma: '', f: 'cos(x)', x0: 0, x1: 2, left: D('dirichlet', '0'), right: D('none') },
    exact: (x) => Math.sin(x) },
  { name: "u'' + 4u = 0, başlangıç değer", inputs: { alpha: '1', beta: '', gamma: '4', f: '0', x0: 0, x1: 1, left: D('both', '1', '0'), right: D('none') },
    exact: (x) => Math.cos(2 * x) }
];
const MESHES = [[1, 1], [2, 1], [3, 2], [5, 3], [7, 4], [16, 2], [64, 4]];

function build(pr, n, p, withReference) {
  const problem = F.buildProblem(pr.inputs), mesh = F.makeMesh(pr.inputs.x0, pr.inputs.x1, n, p);
  const c = F.METHODS.galerkin.solve(problem, mesh).c;   // sayfadaki tek yöntem; başlangıç değer problemini de o çözer
  const reference = withReference && pr.exact ? { u: pr.exact, type: 'analytical', level: 'verified symbolically' } : null;
  return { problem, mesh, c, data: SIM.buildSimulationData({ problem, mesh, c, xs: F.sampleXs(mesh), reference, maxError: reference ? { x: 0.3, value: 1e-3 } : null }) };
}

// 1–2. Veri sözleşmesi = çözücünün verisi; eksik veri null
let coordOk = true, connOk = true, femOk = true, sumWorst = 0, basisWorst = 0, nullOk = true, timeOk = true, coefOk = true;
for (const pr of PROBLEMS) for (const [n, p] of MESHES) {
  let built;
  try { built = build(pr, n, p, true); } catch (err) { check(`${pr.name}, n=${n}, p=${p}: çözülebildi`, false, err.message); continue; }
  const { mesh, c, data } = built, d = data.discretization;
  coordOk &&= d.coordinates.length === mesh.dofs && d.coordinates.every((x, i) => x === mesh.nodeX(i)) && d.nodeCount === n * p + 1 && d.elementCount === n;
  connOk &&= d.connectivity.every((nodes, e) => nodes.length === p + 1 && nodes.every((j, k) => j === e * p + k + 1));
  femOk &&= data.solution.samples.fem.every((v, i) => v === F.evaluate(mesh, c, data.solution.samples.x[i]));
  coefOk &&= data.solution.nodalCoefficients.every((v, i) => v === c[i]) && data.solution.nodalCoefficients !== c;
  // Simülasyonun son yapım adımındaki tam toplam Σ a_j N_j (eleman örneklerinde) = çözücünün u_h'si
  for (let e = 1; e <= n; e++) data.basis.elementX[e - 1].forEach((x, s) => {
    const sum = d.connectivity[e - 1].reduce((acc, j, k) => acc + c[j - 1] * data.basis.N[k][s], 0);
    sumWorst = Math.max(sumWorst, Math.abs(sum - F.evaluate(mesh, c, x)) / Math.max(1, Math.abs(sum)));   // u_h sürekli: sınırda iki eleman aynı değeri verir
  });
  // Destek ve N_j değerleri basis_core.globalBasis ile aynı (küçük ağlarda)
  if (n <= 7) {
    const grid = BASIS.makeGrid(mesh.b - mesh.a, n, p, mesh.a);
    for (let j = 1; j <= d.nodeCount; j++) {
      const g = BASIS.globalBasis(grid, j, 33), sup = data.basis.support[j - 1];
      if (g.support.join() !== sup.map((s) => s.element).join()) basisWorst = Infinity;
      for (const { element, local } of sup) {
        const seg = g.segs[element - 1];
        seg.N.forEach((v, s) => { basisWorst = Math.max(basisWorst, Math.abs(v - data.basis.N[local - 1][s])); });
      }
    }
  }
  const noRef = build(pr, n, p, false).data;
  nullOk &&= noRef.solution.samples.reference === null && noRef.solution.samples.error === null && noRef.metadata.referenceType === null
    && noRef.metadata.errorAvailable === false && noRef.metadata.errorKind === null && noRef.metadata.maxError === null;
  timeOk &&= data.temporal === null && data.metadata.physicalTimeAvailable === false && data.problem.equationOrderT === 0 && data.problem.timeDomain === null;
}
check('düğüm koordinatları = çözücünün mesh.nodeX(i); n_n = n·p + 1', coordOk);
check('bağlantı: e elemanının k. yerel düğümü = (e − 1)·p + k (1 tabanlı)', connOk);
check('u_h örnekleri = çözücünün evaluate(mesh, c, x) (birebir)', femOk);
check('a_j = çözümün c vektörü (birebir, kopya)', coefOk);
check('tam toplam Σ a_j N_j = çözücünün u_h(x) (5 problem × 7 ağ)', sumWorst < 1e-13, sumWorst.toExponential(1));
check('N_j destek ve değerleri = basis_core.globalBasis', basisWorst < 1e-15, String(basisWorst));
check('referans yoksa referans/hata/maxError null (0 ile doldurulmaz)', nullOk);
check('durağan problem: temporal null, fiziksel zaman yok, zaman mertebesi 0', timeOk);

// Analitik ↔ ince ağ etiketi: hata türü doğru adlandırılır
{
  const pr = PROBLEMS[0], problem = F.buildProblem(pr.inputs), mesh = F.makeMesh(0, 1, 4, 2), c = F.METHODS.galerkin.solve(problem, mesh).c;
  const fine = SIM.buildSimulationData({ problem, mesh, c, xs: F.sampleXs(mesh), reference: { u: pr.exact, type: 'fine_mesh', level: 'n = 128, p = 4' } });
  const exact = SIM.buildSimulationData({ problem, mesh, c, xs: F.sampleXs(mesh), reference: { u: pr.exact, type: 'analytical', level: 'verified symbolically' } });
  const sf = SIM.buildStages(fine), se = SIM.buildStages(exact);
  check('ince ağ referansı: "difference", "not the exact error"', fine.metadata.errorKind === 'difference' && sf[7].title === 'Difference from the reference' && /not the exact error/.test(sf[7].text));
  check('analitik referans: "Pointwise error", e = u − u_h', exact.metadata.errorKind === 'exact' && se[7].title === 'Error measures' && /e\(x\) = u\(x\) − u_h\(x\)/.test(se[7].text));
  const holes = SIM.buildSimulationData({ problem, mesh, c, xs: F.sampleXs(mesh), reference: { u: (x) => (x > 0.5 ? NaN : x), type: 'analytical' } });
  const svg = SIM.renderFrame(holes, SIM.buildStages(holes), 1e9, { width: 600 }).svg;
  check('referans tanımsız noktada null (0 değil) ve çizim kırılır', holes.solution.samples.reference.some((v) => v === null)
    && holes.solution.samples.error.every((v, i) => (holes.solution.samples.reference[i] === null) === (v === null)) && !/NaN/.test(svg));
  const noRef = SIM.buildStages(SIM.buildSimulationData({ problem, mesh, c, xs: F.sampleXs(mesh) }));
  check('referans yoksa karşılaştırma ve hata aşamaları hiç yok', noRef.length === 6 && !noRef.some((s) => s.id === 'reference' || s.id === 'error'));
  let threw = 0;
  for (const bad of [c.slice(1), c.map((v, i) => (i === 2 ? NaN : v)), null])
    try { SIM.buildSimulationData({ problem, mesh, c: bad, xs: F.sampleXs(mesh) }); } catch { threw++; }
  check('eksik/yanlış katsayı vektörü reddedilir (uydurma yok)', threw === 3);
}

// 4–6. Her problem × ağ × genişlik × her kare: geçerli SVG, etiket çakışması yok, "time" iddiası yok
let frameOk = true, overlapOk = true, timeWord = true, locateOk = true, partialOk = true, badFrame = '';
for (const pr of PROBLEMS) for (const [n, p] of MESHES) for (const width of [260, 380, 900]) {
  const { data } = build(pr, n, p, true), stages = SIM.buildStages(data);
  const total = SIM.locateStep(stages, Infinity).total;
  locateOk &&= total === stages.reduce((s, st) => s + st.steps, 0) && SIM.locateStep(stages, -5).step === 0 && SIM.locateStep(stages, NaN).step === 0;
  for (let k = 0; k < total; k++) {
    const where = SIM.locateStep(stages, k);
    locateOk &&= where.step === k && stages[where.index].offset + where.local === k;
    const f = SIM.renderFrame(data, stages, k, { width, node: 2, element: n, probe: (data.problem.domain[0] + data.problem.domain[1]) / 2 });
    if (/NaN|undefined|Infinity|null/.test(f.svg) || !f.svg.startsWith('<svg') || !f.svg.endsWith('</svg>')) { frameOk = false; badFrame ||= `${pr.name} n=${n} p=${p} w=${width} k=${k}`; }
    if (/\btime\b|t = /i.test(f.note + stages[where.index].text + stages[where.index].title)) timeWord = false;
    // Etiket çakışması: aynı satırdaki (aynı y) düğüm/eleman numaraları en az kendi genişliği kadar ayrık
    for (const cls of ['sim-node', 'sim-el', 'sim-sel']) {
      const rows = {};
      for (const m of f.svg.matchAll(new RegExp(`<text x="([\\d.]+)" y="([\\d.]+)" text-anchor="middle" class="${cls}">([^<]*)<`, 'g')))
        (rows[m[2]] ||= []).push({ x: Number(m[1]), w: m[3].length * 6.4 });
      for (const row of Object.values(rows)) {
        row.sort((a, b) => a.x - b.x);
        for (let i = 1; i < row.length; i++) if (row[i].x - row[i - 1].x < (row[i].w + row[i - 1].w) / 2) overlapOk = false;
      }
    }
    if (where.stage.id === 'approx' && where.local < where.stage.steps - 1) partialOk &&= /Partial sum of \d+ of \d+ terms/.test(f.note) && /not yet u_h/.test(f.note);
  }
}
check('bütün kareler geçerli SVG (NaN/undefined/Infinity/null yok): 5 problem × 7 ağ × 3 genişlik', frameOk, badFrame);
check('adım ↔ aşama eşlemesi tutarlı; aralık dışı adım kırpılır', locateOk);
check('aynı satırdaki numara etiketleri üst üste binmez (260–900 px)', overlapOk);
check('durağan yapım adımları hiçbir yerde "time" diye etiketlenmez', timeWord);
check('kısmi toplam "u_h" diye etiketlenmez (yalnız son adımda u_h)', partialOk);

// Uç koşulu etiketleri probleme göre doğru (verilmeyen uç: 2. mertebe BVP doğal, BVP dışı serbest)
{
  const lab = (pr, side) => { const d = build(pr, 2, 1, false).data; return SIM.conditionLabel(d.problem.boundaryConditions[side === 'left' ? 0 : 1], d.problem); };
  const natural = { name: 'natural', inputs: { alpha: '1', beta: '', gamma: '', f: '1', x0: 0, x1: 1, left: D('dirichlet', '0'), right: D('none') } };
  check('2. mertebe sınır değer, verilmeyen uç: "natural: α u′ = 0"', lab(natural, 'right') === 'natural: α u′ = 0');
  check('başlangıç değer problemi: "u = 1, u′ = 0" ve serbest uç', lab(PROBLEMS[4], 'left') === 'u = 1, u′ = 0' && lab(PROBLEMS[4], 'right') === 'free end');
  check('1. mertebe, verilmeyen uç serbest; akı koşulu "α u′ = 2"', lab(PROBLEMS[3], 'right') === 'free end' && lab(PROBLEMS[2], 'right') === 'α u′ = 2');
}

// Sorgu: en yakın ÖRNEK noktanın hesaplanmış değerleri
{
  const { data, mesh, c } = build(PROBLEMS[0], 4, 2, true), q = SIM.probeSample(data, 0.4321);
  const xs = data.solution.samples.x, nearest = xs.reduce((b, x, i) => (Math.abs(x - 0.4321) < Math.abs(xs[b] - 0.4321) ? i : b), 0);
  check('sorgu: en yakın örnek, u_h = evaluate, e = u − u_h', q.index === nearest && q.fem === F.evaluate(mesh, c, xs[nearest])
    && q.reference === PROBLEMS[0].exact(xs[nearest]) && q.error === q.reference - q.fem);
  check('sorgu: geçersiz x → null', SIM.probeSample(data, NaN) === null);
}

// 7. Ağır durum süresi (n = 64, p = 4): veri kurma ve en pahalı kare
{
  const t0 = process.hrtime.bigint(), { data } = build(PROBLEMS[0], 64, 4, true), t1 = process.hrtime.bigint();
  const stages = SIM.buildStages(data), total = SIM.locateStep(stages, Infinity).total;
  let worst = 0;
  for (let k = 0; k < total; k++) { const a = process.hrtime.bigint(); SIM.renderFrame(data, stages, k, { width: 900 }); worst = Math.max(worst, Number(process.hrtime.bigint() - a) / 1e6); }
  const buildMs = Number(t1 - t0) / 1e6;
  check('n = 64, p = 4: veri kurma < 300 ms, en pahalı kare < 40 ms', buildMs < 300 && worst < 40, `kurma ${buildMs.toFixed(1)} ms (çözüm dahil), kare ${worst.toFixed(1)} ms`);
}

// 8. Sayfa entegrasyonu: panel ilk sonuç paneli; modül yoksa gizli; CORE'a ve kâğıda bağlanmaz
{
  const main = html.slice(html.indexOf('<main>'), html.indexOf('</main>'));
  const ids = [...main.matchAll(/<section[^>]*id="([^"]+)"/g)].map((m) => m[1]).filter((id) => id !== 'empty-state');
  check('FEM Simulation: sonuç alanının ilk iki paneli Approximation, sonra Solution', ids[0] === 'sim-panel' && ids[1] === 'sol-panel', ids.slice(0, 3).join(' → '));
  check('modül yoksa panel gizli (SIM koruması), betik etiketi var', /const SIM = BASIS && typeof femSimulation === 'function' \?/.test(html)
    && /const visible = !!SIM && !!view;/.test(html) && /\$\(`\$\{pk\}-panel`\)\.hidden = !visible/.test(html) && html.includes('<script src="fem_sim.js"></script>'));
  check('CORE bloğu fem_sim.js\'e başvurmaz', !/femSimulation|fem_sim/.test(src));
  const simSrc = fs.readFileSync(__dirname + '/fem_sim.js', 'utf8');
  check('fem_sim.js DOM\'a dokunmaz ve çözücü çağırmaz', !/document\.|window\.|METHODS|\.solve\(|assemble/.test(simSrc));
  check('zaman çubuğu "not physical time" diye etiketli', /aria-label="Construction step \(not physical time\)"/.test(html) && /not physical time<\/span>/.test(html));
}

// 9. 2B ↔ 3B ortak kare durumu ve 3B modülün sınırları (WebGL node'da yok: davranış tarayıcıda, sınırlar burada sınanır)
{
  let same = true;
  for (const pr of PROBLEMS.slice(0, 2)) for (const [n, p] of [[1, 1], [5, 3], [64, 4]]) {
    const { data } = build(pr, n, p, true), stages = SIM.buildStages(data), total = SIM.locateStep(stages, Infinity).total;
    for (let k = 0; k < total; k++) {
      const fs = SIM.frameState(data, stages, k, { node: 2, element: 1 }), f = SIM.renderFrame(data, stages, k, { width: 600, node: 2, element: 1 });
      if (fs.stage.id === 'approx' && !fs.complete && !f.note.startsWith(`Partial sum of ${fs.terms} of ${data.discretization.nodeCount} terms`)) same = false;
      if (fs.stage.id === 'mesh' && !fs.finalMesh && !f.note.startsWith(`Elements 1–${fs.shownElements} of ${n}`)) same = false;
      if (fs.stage.id === 'global' && fs.local === 2 && !f.note.includes(`nodes ${fs.neighbors.join(', ')}`)) same = false;
    }
  }
  check('2B notu ile ortak kare durumu (frameState) aynı terim/eleman/komşu sayılarını verir', same);
  const r = SIM.valueRange(build(PROBLEMS[0], 4, 2, true).data, 'approx');
  check('valueRange: 0 içeride, sonlu; fonksiyon çizmeyen aşamada null', r.lo < 0 && r.hi > 0 && Number.isFinite(r.lo + r.hi) && SIM.valueRange(build(PROBLEMS[0], 4, 2, true).data, 'mesh') === null);
  const s3 = fs.readFileSync(__dirname + '/fem_sim3d.js', 'utf8'), simSrc = fs.readFileSync(__dirname + '/fem_sim.js', 'utf8');
  const imports = [...s3.matchAll(/^import .* from '([^']+)';$/gm)].map((m) => m[1]);
  check('3B modül yalnız three ve three/addons içe aktarır (çözücü yok)', imports.length && imports.every((m) => m === 'three' || m.startsWith('three/addons/'))
    && !/METHODS|\.solve\(|evaluate\(|buildProblem/.test(s3));
  check('3B modül 2B ile ortak durumu okur (frameState, valueRange), kendi ilerleme kuralı yok', /sim\.frameState\(/.test(s3) && /sim\.valueRange\(/.test(s3)
    && !/Math\.ceil\(\(k \+ 1\)/.test(s3) && /frameState\(data, stages, step, opts\)/.test(simSrc.split('function renderFrame')[1]));
  check('3B modül GPU kaynaklarını bırakır ve sürekli çizim döngüsü kurmaz', /geometry\?\.dispose\(\)/.test(s3) && /m\.dispose\(\)/.test(s3) && /map\?\.dispose\(\)/.test(s3)
    && /renderer\.dispose\(\)/.test(s3) && /observer\.disconnect\(\)/.test(s3) && /controls\.dispose\(\)/.test(s3) && !/setAnimationLoop/.test(s3));
  check('3B derinlik ekseni "not a physical dimension" diye etiketli; hareket azaltma tercihine uyulur', /not a physical dimension/.test(s3) && /prefers-reduced-motion/.test(s3));
  const map = html.match(/<script type="importmap">([^<]+)<\/script>/);
  const imp = map ? JSON.parse(map[1]).imports : {};
  check('Three.js sürümü sabit (0.160.0) ve yalnız import map ile', imp.three === 'https://cdn.jsdelivr.net/npm/three@0.160.0/build/three.module.js'
    && imp['three/addons/'] === 'https://cdn.jsdelivr.net/npm/three@0.160.0/examples/jsm/');
  check('3B yalnız istenince yüklenir (import()); yüklenemezse 2B ve neden yazılır', /await import\('\.\/fem_sim3d\.js'\)/.test(html) && !/<script[^>]+fem_sim3d/.test(html)
    && /3D view unavailable \(\$\{err\.message\}\); showing the 2D view\./.test(html));
  const server = fs.readFileSync(__dirname + '/server.py', 'utf8');
  check('sunucu: 3B modül sabit beyaz listede', /self\.path in \('\/basis_core\.js', '\/fem_sim\.js', '\/fem_sim3d\.js'\)/.test(server));
}

// 10. Girdi fonksiyonları ve akış
{
  const { data, problem } = build(PROBLEMS[2], 4, 2, true);   // α = 1 + x², β yok, γ = exp(x), f = 1
  const inp = data.problem.inputFunctions, xs = data.solution.samples.x;
  check('girdi fonksiyonları: sıfır terim (β) yok, f her zaman var; sırası α, γ, f', inp.map((f) => f.key).join() === 'alpha,gamma,f');
  check('girdi değerleri = çözücünün ayrıştırdığı ifade (problem.coef) birebir', inp.every((fn) => fn.values.every((v, i) => v === problem.coef[fn.key].f(xs[i]))));
  check('sabit girdi işaretli (f = 1), değişken olan null', inp.find((f) => f.key === 'f').constant === 1 && inp.find((f) => f.key === 'alpha').constant === null);
  const zeroF = build({ name: 'f0', inputs: { alpha: '1', beta: '', gamma: '', f: '', x0: 0, x1: 1, left: D('dirichlet', '0'), right: D('dirichlet', '1') } }, 2, 1, false).data;
  check('f = 0 da çizilir (bilgidir), α sabit 1', zeroF.problem.inputFunctions.map((f) => `${f.key}=${f.constant}`).join() === 'alpha=1,f=0');
  const stages = SIM.buildStages(data), f2 = SIM.renderFrame(data, stages, 0, { width: 700, probe: 0.5 });
  check('1. aşama 3 adım; 1. adımda (panel açılışında) her girdi için bir panel ve değer okuma', stages[0].steps === 3 && (f2.svg.match(/class="sim-fn"/g) || []).length === inp.length
    && f2.scale.ranges.length === inp.length && f2.probe.inputs.length === inp.length && f2.probe.inputs.every((q) => Number.isFinite(q.value)));
  const narrow = SIM.renderFrame(data, stages, 0, { width: 300 });
  check('dar ekranda paneller 2 sütun, sayfa dışına taşmaz', narrow.scale.ranges.every(([l, r]) => l >= 0 && r <= 300) && new Set(narrow.scale.ranges.map(([l]) => Math.round(l))).size <= 2);
  let flowOk = true, ghostOk = true, which = '';
  for (const pr of PROBLEMS) for (const [n, p] of [[1, 1], [4, 2], [9, 3]]) {
    const d = build(pr, n, p, true).data, st = SIM.buildStages(d), total = SIM.locateStep(st, Infinity).total;
    for (let k = 0; k < total; k++) {
      const fs = SIM.frameState(d, st, k, {}), f = SIM.renderFrame(d, st, k, { width: 640 });
      if (!/class="(sweep|enter)"/.test(f.svg)) { flowOk = false; which ||= `${pr.name} n=${n} p=${p} k=${k} ${fs.stage.id}`; }
      if (fs.stage.id === 'approx' && !fs.complete) ghostOk &&= (fs.prevTerms > 0) === /<g opacity="0\.28">/.test(f.svg) && fs.prevTerms === fs.newestFrom - 1;
    }
  }
  check('akış: her karede bu adımda yeni gelen bir öğe canlandırılır', flowOk, which);
  check('önceki kısmi toplam izi yalnız önceki adımda terim varsa ve tam o terim sayısıyla', ghostOk);
  const css = html.match(/@keyframes sim-sweep[^}]*}[^}]*}/)?.[0] || '';
  check('akış yalnız görünme biçimi (clip/opacity), hareket azaltma tercihinde kapalı', /clip-path: inset/.test(css)
    && /prefers-reduced-motion: reduce\) \{ \.sim-chart \.sweep, \.sim-chart \.enter \{ animation: none; \}/.test(html));
}

// 11. Radyan ekseni (trigonometrik girdi): işaretlerin yazılışı; çizilen veri ve örnekleme değişmez
{
  const { usesTrigonometry: trig, angularTicks: ticks } = require('./fem_sim.js');
  const labels = (lo, hi, c = 5) => (ticks(lo, hi, c) || []).map((t) => t.label).join(' ');
  check('trigonometri algısı: sin(x), cos(2*x), tan(x) evet', trig(['', '', '', 'sin(x)']) && trig(['cos(2*x)']) && trig(['1 + tan(x/2)']) && trig(['SIN (x)']));
  check('açı olmayanlar hayır: sin(pi*x), sinh, asin, cosh, sinx, boş', !trig(['sin(pi*x)']) && !trig(['sinh(x)']) && !trig(['asin(x)', 'cosh(x)']) && !trig(['2x + x^2', '']) && !trig([undefined, null]) && !trig(['mysin(x)']));
  check('iç içe parantez: sin((pi*x)) hayır, sin(x + pi) açı (π toplanan, çarpan değil sayılmaz → evet)', !trig(['sin((pi*x))']) === true && trig(['sin(x)*cos(pi*x)']));
  check('[0, 2π]: 0 π/2 π 3π/2 2π', labels(0, 2 * Math.PI) === '0 π/2 π 3π/2 2π', labels(0, 2 * Math.PI));
  check('[−π, π]: −π −π/2 0 π/2 π', labels(-Math.PI, Math.PI) === '−π −π/2 0 π/2 π', labels(-Math.PI, Math.PI));
  check('[0, π]: 4 işaret (π/3 adımı sığmaz, π/4 ya da π/2)', /^0 .*π$/.test(labels(0, Math.PI)) && ticks(0, Math.PI, 5).length >= 3, labels(0, Math.PI));
  check('uç değer tam π katıyken kırpılmaz (yuvarlama payı)', ticks(0, 6.283185307179586, 5).at(-1).label === '2π' && ticks(0, 3.141592653589793 * 2 + 1e-12, 5).at(-1).label === '2π');
  check('[0, 1] (π katı yok): null → onluk işaretlere dönülür', ticks(0, 1) === null && ticks(0.2, 1.5) === null);
  check('[0, 10π]: adım büyür (≤ count), etiketler tam sayı π', ticks(0, 10 * Math.PI, 6).length <= 6 && ticks(0, 10 * Math.PI, 6).every((t) => /^(0|\d+π)$/.test(t.label)), labels(0, 10 * Math.PI, 6));
  check('π/3 adımı: 2π/3 ve π doğru kesir', /2π\/3/.test(labels(0, Math.PI, 7)) || /π\/3/.test(labels(0, Math.PI, 7)) || /π\/4/.test(labels(0, Math.PI, 7)), labels(0, Math.PI, 7));
  check('geçersiz aralık: null (NaN, tersi, eşit, sonsuz)', ticks(NaN, 1) === null && ticks(2, 1) === null && ticks(1, 1) === null && ticks(0, Infinity) === null);
  let ordered = true, valuesOk = true;
  for (const [lo, hi] of [[0, 7], [-10, 10], [0.5, 40], [-1, 6.5], [100, 120], [-1e3, 1e3]]) for (const c of [3, 4, 6, 8]) {
    const t = ticks(lo, hi, c);
    if (!t) continue;
    ordered &&= t.length <= c && t.every((q, i) => i === 0 || q.value > t[i - 1].value) && t.every((q) => q.value >= lo - 1e-9 && q.value <= hi + 1e-9);
    valuesOk &&= t.every((q) => Math.abs(q.value / Math.PI - Math.round(q.value / Math.PI * 12) / 12) < 1e-9 || Math.abs(q.value / Math.PI - Math.round(q.value / Math.PI)) < 1e-9);
  }
  check('çeşitli aralık ve sayıda: sıralı, aralık içinde, ≤ count, değerler π/12 katı', ordered && valuesOk);
  const prob = build(PROBLEMS[0], 4, 2, true), mk = (angular) => SIM.buildSimulationData({ problem: prob.problem, mesh: prob.mesh, c: prob.c, xs: F.sampleXs(prob.mesh), angular });
  const on = mk(true), off = mk(false);
  check('radyan bayrağı veriyi değiştirmez: örnekler, u_h ve katsayılar aynı', on.problem.angularX === true && off.problem.angularX === false
    && JSON.stringify(on.solution) === JSON.stringify(off.solution) && JSON.stringify(on.discretization) === JSON.stringify(off.discretization));
  const wide = { ...PROBLEMS[0], inputs: { ...PROBLEMS[0].inputs, x0: 0, x1: 2 * Math.PI } }, wp = F.buildProblem(wide.inputs), wm = F.makeMesh(0, 2 * Math.PI, 4, 2), wc = F.METHODS.galerkin.solve(wp, wm).c;
  const wd = SIM.buildSimulationData({ problem: wp, mesh: wm, c: wc, xs: F.sampleXs(wm), angular: true }), ws = SIM.buildStages(wd);
  const frame = SIM.renderFrame(wd, ws, SIM.locateStep(ws, Infinity).total - 1, { width: 700 }), early = SIM.renderFrame(wd, ws, 4, { width: 700 });
  check('2B: radyanda π etiketleri ve "rad →"; onlukta π yok', />2π\/3</.test(early.svg) && />rad →</.test(early.svg) && !/π\/2</.test(SIM.renderFrame(off && SIM.buildSimulationData({ problem: wp, mesh: wm, c: wc, xs: F.sampleXs(wm), angular: false }), ws, 4, { width: 700 }).svg));
  check('3B ve 2B aynı işaretleri kullanır (sim.xTicks)', SIM.xTicks(0, 2 * Math.PI, 5, true).map((t) => t.label).join(' ') === '0 π/2 π 3π/2 2π' && SIM.xTicks(0, 1, 5, true)[0].label === '0' && SIM.xTicks(0, 1, 5, false).length >= 3);
  check('sayfa: radyan yalnız eksen yazılışı; yakınlaştırma/log eksende kapalı, işleyici koruması var', /const angular = cfg\.xAngular && !cfg\.xLog && !cfg\.xTicks && typeof angularTicks === 'function'/.test(html)
    && /xAngular: radianAxis\(\)/.test(html) && !/xAngular: radianAxis\(\)[^\n]*conv/.test(html));
}

// 12. Galerkin cebiri paneli (Solution): veri çözücüyle ve sunucuyla tutar mı; aşamalar, formüller, kenar durumlar
{
  const MESHES2 = [[1, 1], [2, 2], [4, 2], [3, 3], [5, 4], [16, 2]];
  const PR2 = [...PROBLEMS, { name: 'doğal koşul + taşınım', inputs: { alpha: '1', beta: '0.5', gamma: '', f: '1', x0: 0, x1: 1, left: D('none'), right: D('dirichlet', '0') } }];
  const buildA = (pr, n, p, withRef = true) => {
    const problem = F.buildProblem(pr.inputs), mesh = F.makeMesh(pr.inputs.x0, pr.inputs.x1, n, p), res = F.METHODS.galerkin.solve(problem, mesh);
    const reference = withRef && pr.exact ? { u: pr.exact, type: 'analytical', level: 'verified symbolically' } : null;
    const data = SIMA.buildSimulationData({ problem, mesh, c: res.c, xs: F.sampleXs(mesh), reference, solved: { K: res.K, F: res.F, condition: res.diagnostics.condition },
      tex: { equation: 'u\'\' = f', conditions: 'u(0) = 0' } });
    return { problem, mesh, res, data };
  };
  let dK = 0, dF = 0, dRed = 0, dArea = 0, dim = true, frames = 0, badFrames = 0, badFormula = 0, whichBad = '';
  for (const pr of PR2) for (const [n, p] of MESHES2) {
    let b;
    try { b = buildA(pr, n, p); } catch (err) { continue; }
    const A = b.data.algebra, N = b.mesh.dofs, Ks = Array.from({ length: N }, () => new Array(N).fill(0)), Fs = new Array(N).fill(0);
    for (let e = 0; e < n; e++) for (let i = 0; i <= p; i++) { Fs[e * p + i] += A.elementF[e][i]; for (let j = 0; j <= p; j++) Ks[e * p + i][e * p + j] += A.elementK[e][i][j]; }
    const scale = Math.max(1, ...b.res.K.flat().map(Math.abs));
    for (let i = 0; i < N; i++) { for (let j = 0; j < N; j++) dK = Math.max(dK, Math.abs(Ks[i][j] - b.res.K[i][j]) / scale); dF = Math.max(dF, Math.abs(Fs[i] + A.Fb[i] - b.res.F[i]) / Math.max(1, Math.abs(b.res.F[i]))); }
    const sol = F.solveLinear(A.reducedK, A.reducedF, {});
    A.free.forEach((idx, q) => { dRed = Math.max(dRed, Math.abs(sol[q] - b.res.c[idx]) / Math.max(1, Math.abs(b.res.c[idx]))); });
    dim &&= A.reducedK.length === A.rows.length && A.reducedK.every((r) => r.length === A.free.length) && A.reducedF.length === A.rows.length && A.rows.length === A.free.length;
    if (n <= 5) for (let e = 0; e < n; e++) for (let i = 0; i <= p; i++) for (let j = 0; j <= p; j++) {   // integrand alanı = K⁽ᵉ⁾_ij (Simpson, 32 aralık)
      const xs = b.data.basis.elementX[e], hh = (xs[xs.length - 1] - xs[0]) / (xs.length - 1);
      const tot = (q) => [0, 1, 2].reduce((a, t) => a + A.terms[((((e * (p + 1) + i) * (p + 1) + j) * 3 + t) * A.termSamples) + q], 0);
      let area = 0;
      for (let q = 0; q < xs.length; q++) area += (q === 0 || q === xs.length - 1 ? 1 : q % 2 ? 4 : 2) * tot(q);
      area *= hh / 3;
      dArea = Math.max(dArea, Math.abs(area - A.elementK[e][i][j]) / Math.max(1, Math.abs(A.elementK[e][i][j])));
    }
    const stages = SIMA.buildStages(b.data);
    for (const panel of ['approx', 'solution']) {
      const sub = SIMA.stagesFor(stages, panel), total = SIMA.locateStep(sub, Infinity).total;
      for (const width of [300, 700]) for (let k = 0; k < total; k++) {
        const f = SIMA.renderFrame(b.data, sub, k, { width, cell: [1, 1], node: 2, element: n }); frames++;
        if (/NaN|undefined|Infinity|null/.test(f.svg)) { badFrames++; whichBad ||= `${pr.name} n=${n} p=${p} ${panel} k=${k} w=${width}`; }
        if (!Array.isArray(f.formulas) || f.formulas.some((t) => /NaN|undefined|null/.test(t) || !t.trim())) badFormula++;
      }
    }
  }
  check('Σ K⁽ᵉ⁾ = çözücünün K matrisi (7 problem × 6 ağ, işaret kuralı dahil)', dK < 1e-13, dK.toExponential(1));
  check('Σ F⁽ᵉ⁾ + F_b = çözücünün F vektörü', dF < 1e-13, dF.toExponential(1));
  check('indirgenmiş sistemin çözümü = çözücünün a katsayıları (başlangıç değer ve 1. mertebe dahil)', dRed < 1e-12 && dim, dRed.toExponential(1));
  check('integrand eğrisinin altındaki alan = K⁽ᵉ⁾_ij (Simpson, 7 problem × 6 ağ; çizim çekirdek kuadratürüyle tutar)', dArea < 2e-3, dArea.toExponential(1));
  check('her iki panelin her karesi geçerli SVG ve formüllü', badFrames === 0 && badFormula === 0, `${frames} kare ${whichBad}`);

  const b0 = buildA(PROBLEMS[0], 4, 2), st0 = SIMA.buildStages(b0.data);
  const ids = (panel) => SIMA.stagesFor(st0, panel).map((q) => q.id).join(' ');
  check('Approximation paneli: problem, ağ, eleman, yerel şekil, global taban', ids('approx') === 'domain mesh element local global');
  check('Solution paneli = kâğıdın bölümleri sırayla: yaklaşım → residual → zayıf form → eleman → montaj → sınır → çözme → yeniden kurma → karşılaştırma → hata',
    ids('solution') === 'ansatz residual weak elementK assembly bc solve approx reference error', ids('solution'));
  const sub = SIMA.stagesFor(st0, 'solution');
  check('paneller kendi içinde 0\'dan numaralanır ve toplam adım = aşama adımları toplamı', SIMA.stagesFor(st0, 'approx')[0].offset === 0 && sub[0].offset === 0
    && SIMA.locateStep(sub, Infinity).total === sub.reduce((a, q) => a + q.steps, 0));
  const o1 = SIMA.buildStages(buildA(PROBLEMS[3], 3, 2).data), o2 = st0;
  check('2. mertebede zayıf form 2 adım, bc 4 adım; 1. mertebede zayıf form aşaması yok (kısmi integrasyon yok), bc 3 adım',
    o2.find((q) => q.id === 'weak').steps === 2 && o2.find((q) => q.id === 'bc').steps === 4 && !o1.some((q) => q.id === 'weak') && o1.find((q) => q.id === 'bc').steps === 3);
  check('elementSystem yoksa cebir aşamaları hiç yok (uydurma matris yok), yeniden kurma ve hata kalır',
    SIM.buildStages(SIM.buildSimulationData({ problem: b0.problem, mesh: b0.mesh, c: b0.res.c, xs: F.sampleXs(b0.mesh), solved: { K: b0.res.K, F: b0.res.F },
      reference: { u: PROBLEMS[0].exact, type: 'analytical' } })).map((q) => q.id).join(' ') === 'domain mesh element local global approx reference error');

  // formüller: yalnız sıfır olmayan terimler; sabit kalıplar kâğıtla aynı
  const eqs = (pr, k, id, opts = {}) => { const bb = buildA(pr, 3, 2), stg = SIMA.stagesFor(SIMA.buildStages(bb.data), 'solution'), ix = stg.findIndex((q) => q.id === id);
    return SIMA.stageFormulas(bb.data, SIMA.frameState(bb.data, stg, stg[ix].offset + k, opts), opts).join(' '); };
  const f0 = eqs(PROBLEMS[0], 0, 'elementK'), f1 = eqs(PROBLEMS[3], 0, 'elementK');
  check('eleman formülü yalnız denklemdeki terimleri içerir (β yok → β terimi yok; α yok → α terimi yok)', /\\alpha/.test(f0) && /\\gamma/.test(f0) && !/\\beta/.test(f0) && /\\beta/.test(f1) && !/\\alpha/.test(f1) && !/\\gamma/.test(f1));
  check('2. mertebede F⁽ᵉ⁾_i = −∫…, 1. mertebede +∫… (çekirdeğin işaret kuralı); zayıf formda F_b yalnız 2. mertebede', /F\^\{\(e\)\}_i = -\\int/.test(eqs(PROBLEMS[0], 0, 'elementK')) && /F\^\{\(e\)\}_i = \\int/.test(eqs(PROBLEMS[3], 0, 'elementK')) && !/F\^\{\(e\)\}_i = -/.test(eqs(PROBLEMS[3], 0, 'elementK')) && /F_\{b,i\}/.test(eqs(PROBLEMS[0], 1, 'weak')));
  check('seçili hücrenin sayısı K⁽ᵉ⁾_ij formülde çekirdek değeriyle yazılır', (() => { const bb = buildA(PROBLEMS[0], 3, 2), v = bb.data.algebra.elementK[1][0][1]; return new RegExp(`K\\^\\{\\(2\\)\\}_\\{12\\} = ${String(Number(v.toPrecision(4))).replace('-', '-').replace('.', '\\.')}`).test(eqs(PROBLEMS[0], 2, 'elementK', { element: 2 })); })());

  // kâğıt tek kaynak: sabit kalıplar kâğıt çıktısında aynen bulunur; eleman matrisleri sunucunun sembolik/tam hesabıyla tutar
  const sheetOf = (pr, n, p) => { const r = spawnSync('python3', ['fem_sheet.py'], { cwd: __dirname, input: JSON.stringify({ ...pr.inputs, n, degree: p, symbolic: false }), encoding: 'utf8', timeout: 120000 }); return JSON.parse(r.stdout); };
  const sheet = sheetOf(PROBLEMS[0], 4, 2), jsA = buildA(PROBLEMS[0], 4, 2).data;
  const solStages = SIMA.stagesFor(SIMA.buildStages(jsA), 'solution'), formulasOf = (id, k) => SIMA.stageFormulas(jsA, SIMA.frameState(jsA, solStages, solStages.find((q) => q.id === id).offset + k, {}), {}).join(' ');
  const jsText = [formulasOf('ansatz', 0), formulasOf('assembly', 0), formulasOf('bc', 2)].join(' ');
  for (const frag of ['u_h(x) = \\sum_{j=1}^{9} a_j\\,N_j(x)', 'K_{IJ} = \\sum_{e} K^{(e)}_{ij}', 'F_I = \\sum_e F^{(e)}_i', 'a_{1} = 0\\ \\text{(essential)}'])
    check(`kâğıttaki formül JS altyazısıyla aynı dizge: ${frag.slice(0, 34)}…`, sheet.text.includes(frag) && jsText.includes(frag));
  let dSrv = 0;
  for (const [pr, n, p] of [[PROBLEMS[0], 4, 2], [PROBLEMS[2], 3, 3], [PROBLEMS[3], 3, 2], [PROBLEMS[1], 2, 1]]) {
    const sh = sheetOf(pr, n, p), js = buildA(pr, n, p).data.algebra;
    sh.element_K.forEach((Ke, e) => Ke.forEach((row, i) => row.forEach((v, j) => { dSrv = Math.max(dSrv, Math.abs(v - js.elementK[e][i][j]) / Math.max(1, Math.abs(v))); })));
    sh.element_F.forEach((Fe, e) => Fe.forEach((v, i) => { dSrv = Math.max(dSrv, Math.abs(v - js.elementF[e][i]) / Math.max(1, Math.abs(v))); }));
  }
  check('eleman matrisleri: JS (kuadratür) = sunucu (SymPy, tam) — iki bağımsız çekirdek, 4 problem', dSrv < 1e-9, dSrv.toExponential(1));

  // sayfa: iki panel, panel başına zaman çubuğu, matris aşamaları 3B'ye gönderilmez
  check('sayfa: iki panel, her birinin kendi durumu ve zamanlayıcısı', /id="sim-panel"/.test(html) && /id="sol-panel"/.test(html) && /panels: \{ sim: \{ step: 0, playing: false, speed: 1, view: '2d' \}, sol:/.test(html)
    && /const simTimers = \{ sim: null, sol: null \}/.test(html));
  check('matris aşamaları yalnız 2B: 3B desteklenen aşama listesinde cebir aşaması yok', /const SIM_3D_STAGES = new Set\(\[[^\]]*\]\)/.test(html) && !/SIM_3D_STAGES = new Set\(\[[^\]]*(residual|weak|elementK|assembly|'bc'|solve)/.test(html));
}

// 13. Solution paneli = çözüm kâğıdının adım adım gezintisi: kâğıt blokları, hücre hücre dolan matris, boş K'de sıfır yok
{
  const { parseSheet: parse, sheetPick: pick } = require('./fem_sim.js');
  const WALK = [
    { name: 'sabit katsayı', pr: PROBLEMS[0], n: 4, p: 2 },
    { name: 'değişken katsayı + akı', pr: PROBLEMS[2], n: 3, p: 2 },
    { name: '1. mertebe', pr: PROBLEMS[3], n: 3, p: 2 },
    { name: 'kırpılmış kâğıt (n = 10)', pr: PROBLEMS[0], n: 10, p: 1 }
  ];
  const sheetFor = (pr, n, p, symbolic = true) => { const r = spawnSync('python3', ['fem_sheet.py'], { cwd: __dirname, input: JSON.stringify({ ...pr.inputs, n, degree: p, symbolic }), encoding: 'utf8', timeout: 240000 }); return JSON.parse(r.stdout); };
  const mk = (pr, n, p) => { const problem = F.buildProblem(pr.inputs), mesh = F.makeMesh(pr.inputs.x0, pr.inputs.x1, n, p), res = F.METHODS.galerkin.solve(problem, mesh);
    return { problem, mesh, res, data: SIMA.buildSimulationData({ problem, mesh, c: res.c, xs: F.sampleXs(mesh), solved: { K: res.K, F: res.F, condition: res.diagnostics.condition },
      reference: pr.exact ? { u: pr.exact, type: 'analytical', level: 'verified symbolically' } : null, tex: { equation: 'u\'\' = f', conditions: 'u(0) = 0' } }) }; };
  for (const { name, pr, n, p } of WALK) {
    const sheet = sheetFor(pr, n, p), sec = parse(sheet.text), b = mk(pr, n, p), sub = SIMA.stagesFor(SIMA.buildStages(b.data), 'solution');
    check(`kâğıt ayrıştırma (${name}): bölümler başlıktan tanınır; yeniden kurma, yaklaşım bölümünün üstüne yazmaz`,
      !!sec && ['approx', 'galerkin', 'element', 'assembly', 'bc', 'solution', 'recon'].every((q) => sec[q] && sec[q].blocks.length)
      && sec.approx.title.includes('Approximate solution') && !sec.approx.title.includes('Reconstruction') && sec.recon.title.includes('Reconstruction') && (pr.inputs.alpha === '' ? !sec.weak : !!sec.weak));
    let missing = 0, total = 0, wrong = 0;
    for (const st of sub) for (let k = 0; k < st.steps; k++) for (const el of [1, Math.ceil(n / 2), n]) {
      const f = SIMA.renderFrame(b.data, sub, st.offset + k, { width: 700, sheet: sec, element: el }); total++;
      if (!f.formulaSource || !f.formulas.length) missing++;
      if (f.formulas.some((t) => !t.trim() || /undefined|NaN/.test(t))) wrong++;
    }
    check(`her adımda kâğıdın kendi blokları var (${name}): ${total} adım×eleman, yedek formüle düşen 0`, missing === 0 && wrong === 0, `eksik ${missing}, bozuk ${wrong}`);
  }

  // eleman matrisi adımları: hücre hücre dolar, doğru hücrenin kâğıt bloğu gelir
  {
    const { pr, n, p } = WALK[1], sheet = sheetFor(pr, n, p), sec = parse(sheet.text), b = mk(pr, n, p), sub = SIMA.stagesFor(SIMA.buildStages(b.data), 'solution');
    const st = sub.find((q) => q.id === 'elementK'), P1 = p + 1;
    check('eleman matrisi adım sayısı = 1 + (p+1)² (K hücreleri) + (p+1) (F girdileri) + 1', st.steps === 1 + P1 * P1 + P1 + 1);
    let order = true, counts = true, labels = true;
    for (let k = 0; k < st.steps; k++) {
      const f = SIMA.renderFrame(b.data, sub, st.offset + k, { width: 700, sheet: sec, element: 2 }), texts = (f.svg.match(/class="sim-cell"/g) || []).length;
      const want = k === 0 ? 0 : k <= P1 * P1 ? k : k <= P1 * P1 + P1 ? P1 * P1 + (k - P1 * P1) : P1 * P1 + P1;
      counts &&= texts === want;
      if (k >= 1 && k <= P1 * P1) { const i = Math.floor((k - 1) / P1), j = (k - 1) % P1; labels &&= f.formulas.some((t) => new RegExp(`^\\\\begin\\{aligned\\} K\\^\\{\\((\\d)\\)\\}_\\{${i + 1}${j + 1}\\}`).test(t)); }
      if (k > P1 * P1 && k <= P1 * P1 + P1) { const i = k - P1 * P1 - 1; labels &&= f.formulas.some((t) => new RegExp(`^\\\\begin\\{aligned\\} F\\^\\{\\(\\d\\)\\}_\\{l,${i + 1}\\}`).test(t)); }
    }
    check('matris adım adım dolar: k. adımda k hücre değeri yazılıdır, F girdileri sırayla, son adımda hepsi', counts);
    check('her hücre adımında o hücrenin (i, j) kâğıt bloğu gösterilir (K⁽ᵉ⁾_ij ve F⁽ᵉ⁾_l,i)', labels);
    const f2 = SIMA.renderFrame(b.data, sub, st.offset + 1, { width: 700, sheet: sec, element: 2 });
    check('değişken katsayıda seçili elemanın bloğu seçilir (K⁽²⁾₁₁)', f2.formulas.some((t) => /K\^\{\(2\)\}_\{11\}/.test(t)));
    const cst = WALK[0], sec0 = parse(sheetFor(cst.pr, cst.n, cst.p).text), b0 = mk(cst.pr, cst.n, cst.p), sub0 = SIMA.stagesFor(SIMA.buildStages(b0.data), 'solution'), st0 = sub0.find((q) => q.id === 'elementK');
    const fk = SIMA.renderFrame(b0.data, sub0, st0.offset + 1, { width: 700, sheet: sec0, element: 3 }), ff = SIMA.renderFrame(b0.data, sub0, st0.offset + 1 + 9 + 1, { width: 700, sheet: sec0, element: 3 });
    check('sabit katsayıda K kâğıtta yalnız e = 1 açılır ("tüm elemanlar için aynı"): K bloğu e = 1, yük bloğu seçili elemandan (e = 3)',
      fk.formulas.some((t) => /K\^\{\(1\)\}_\{11\}/.test(t)) && ff.formulas.some((t) => /F\^\{\(3\)\}_\{l,2\}/.test(t)), JSON.stringify(ff.formulas.map((t) => t.slice(0, 40))));
  }

  // boş K'da sıfır yazılmaz; montaj örneği vurgulanır
  {
    const b = mk(PROBLEMS[0], 4, 2), sub = SIMA.stagesFor(SIMA.buildStages(b.data), 'solution'), st = sub.find((q) => q.id === 'assembly');
    const f0 = SIMA.renderFrame(b.data, sub, st.offset, { width: 700 }), f2 = SIMA.renderFrame(b.data, sub, st.offset + 2, { width: 700 });
    check('montajın ilk adımında K boş: hücrelerde "0" yazılı değil, notta "her giriş 0" denir', !/class="sim-cell">[-−]?0</.test(f0.svg) && (f0.svg.match(/class="sim-cell"/g) || []).every((_, i, arr) => arr.length === b.data.discretization.nodeCount) && /every entry 0/.test(f0.note));   // yalnız a₁…a_N simgeleri var
    check('montaj 2. adımda ortak düğüm girişi vurgulanır, değerli hücreler yazılır', /K₍ₚ₊₁,ₚ₊₁₎/.test(f2.note) && f2.svg.includes('stroke:var(--ink)" stroke-width="2"'));
  }

  // sayısal yedek kâğıtta integral adımları yoktur: hücre adımları yedek formüle düşer, diğer bölümler kâğıttan gelir
  {
    const { pr, n, p } = WALK[0], sec = parse(sheetFor(pr, n, p, false).text), b = mk(pr, n, p), sub = SIMA.stagesFor(SIMA.buildStages(b.data), 'solution'), st = sub.find((q) => q.id === 'elementK');
    const cell = SIMA.renderFrame(b.data, sub, st.offset + 1, { width: 700, sheet: sec, element: 2 }), asm = SIMA.renderFrame(b.data, sub, sub.find((q) => q.id === 'assembly').offset, { width: 700, sheet: sec });
    check('sayısal yedek kâğıt: hücre integrali yazılmamıştır → yedek formül (kaynak null), montaj bölümü yine kâğıttan', cell.formulaSource === null && cell.formulas.length > 0 && !!asm.formulaSource);
  }
  const none = SIMA.renderFrame(mk(PROBLEMS[0], 3, 2).data, SIMA.stagesFor(SIMA.buildStages(mk(PROBLEMS[0], 3, 2).data), 'solution'), 0, { width: 700, sheet: null });
  check('kâğıt yoksa (sunucu yok) yerleşik formül, kaynak null', none.formulaSource === null && none.formulas.length > 0);
  check('kâğıt yoksa pick null döner; boş bölüm sözlüğü de null', pick(null, { id: 'bc' }, { local: 0, selEl: 1 }, { order: 2, p: 1, n: 2 }) === null && parse('metin') === null && parse('') === null);
  check('sayfa: eskimiş kâğıt gösterilmez (anahtar tutmuyorsa null), kâğıt yalnız panel görününce istenir, hücre tıklaması adıma bağlı',
    /const currentSimSheet = \(\) => \(simSheet\.status === 'ok' && simSheet\.key === JSON\.stringify\(sheetPayload\(\)\) \? simSheet\.sections : null\)/.test(html)
    && /new IntersectionObserver\(/.test(html) && /simSheetWanted/.test(html) && /ps\.step = st\.offset \+ 1 \+ ci \* \(simCache\.data\.discretization\.polynomialDegree \+ 1\) \+ cj/.test(html));
}

// 14. Panel 1 (Approximation) formülleri de çözüm kâğıdından: Given / Mesh / Shape functions / Approximate solution bölümleri
{
  const { parseSheet: parse } = require('./fem_sim.js');
  const sheetFor = (pr, n, p, symbolic = true) => { const r = spawnSync('python3', ['fem_sheet.py'], { cwd: __dirname, input: JSON.stringify({ ...pr.inputs, n, degree: p, symbolic }), encoding: 'utf8', timeout: 240000 }); return JSON.parse(r.stdout); };
  const mk = (pr, n, p) => { const problem = F.buildProblem(pr.inputs), mesh = F.makeMesh(pr.inputs.x0, pr.inputs.x1, n, p), res = F.METHODS.galerkin.solve(problem, mesh);
    return { data: SIMA.buildSimulationData({ problem, mesh, c: res.c, xs: F.sampleXs(mesh), solved: { K: res.K, F: res.F, condition: res.diagnostics.condition },
      reference: null, tex: { equation: 'u\'\' = f', conditions: 'u(0) = 0' } }) }; };
  const CASES = [
    { name: 'sabit katsayı P2', pr: PROBLEMS[0], n: 4, p: 2 }, { name: 'değişken katsayı + akı P3', pr: PROBLEMS[2], n: 3, p: 3 },
    { name: '1. mertebe P1', pr: PROBLEMS[3], n: 3, p: 1 }, { name: 'başlangıç değer P4', pr: PROBLEMS[4], n: 3, p: 4 },
    { name: 'kırpılmış kâğıt n = 10 P2', pr: PROBLEMS[0], n: 10, p: 2 }
  ];
  const sheets = new Map(CASES.map((c) => [c.name, parse(sheetFor(c.pr, c.n, c.p).text)]));
  for (const c of CASES) {
    const sec = sheets.get(c.name), b = mk(c.pr, c.n, c.p), sub = SIMA.stagesFor(SIMA.buildStages(b.data), 'approx');
    check(`kâğıt bölümleri tanınır (${c.name}): Given, Notation, Mesh, Shape functions, Approximate solution`,
      ['given', 'notation', 'mesh', 'shape', 'approx'].every((q) => sec[q] && sec[q].blocks.length) && sec.shape.title.includes('Shape functions') && sec.given.title === 'Given');
    let missing = 0, bad = 0, total = 0;
    for (const st of sub) for (let k = 0; k < st.steps; k++) for (const el of [1, Math.ceil(c.n / 2), c.n]) {
      const f = SIMA.renderFrame(b.data, sub, st.offset + k, { width: 700, sheet: sec, element: el }); total++;
      if (!f.formulaSource || !f.formulas.length) missing++;
      if (f.formulas.some((t) => !t.trim() || /undefined|NaN/.test(t))) bad++;
    }
    check(`panel 1 her adımda kâğıdın blokları (${c.name}): ${total} adım×eleman, yedeğe düşen 0`, missing === 0 && bad === 0, `eksik ${missing}, bozuk ${bad}`);
  }
  // içerik: doğru adım doğru bloğu getirir
  {
    const c = CASES[0], sec = sheets.get(c.name), b = mk(c.pr, c.n, c.p), sub = SIMA.stagesFor(SIMA.buildStages(b.data), 'approx');
    const at = (id, k, o = {}) => SIMA.renderFrame(b.data, sub, sub.find((q) => q.id === id).offset + k, { width: 700, sheet: sec, ...o });
    const join = (f) => f.formulas.join('\n');
    check('Continuous problem: 1. adım L operatörü, 2. adım denklem + aralık, 3. adım koşullar (kâğıdın Given bölümü)',
      /^\\mathcal\{L\}\(u\)/.test(at('domain', 0).formulas[0]) && /\\in \\Omega = \(/.test(join(at('domain', 1))) && /condition\)/.test(join(at('domain', 2))) && at('domain', 0).formulaSource === 'Given');
    check('Mesh: ilk adım n, h, p (Given); son adım düğüm vektörü ve yerel→global tablosu (Mesh bölümü)',
      /^n = /.test(at('mesh', 0).formulas[0]) && /\\big\[x_j\\big\]/.test(join(at('mesh', sub.find((q) => q.id === 'mesh').steps - 1)))
      && /\\mapsto/.test(join(at('mesh', sub.find((q) => q.id === 'mesh').steps - 1))) && /Mesh/.test(at('mesh', 1).formulaSource));
    check('Element structure: referans eleman s ∈ [0, h]; 2. adım x = x_e + s ve yerel→global tablosu', /Reference element: } s \\in \[0, h\]/.test(join(at('element', 0))) && /^x = x_e \+ s/.test(at('element', 1).formulas[0]) && /\\Omega_e/.test(join(at('element', 1))));
    const L = (k) => join(at('local', k, { element: 1 }));
    check('Local shape functions: k. adımda yalnız N_k(s) satırı (öbürleri yok) ve δ_ij matrisi',
      /N_\{1\}\(s\)/.test(L(0)) && !/N_\{2\}\(s\)/.test(L(0)) && /N_\{2\}\(s\)/.test(L(1)) && !/N_\{1\}\(s\)/.test(L(1)) && /N_\{3\}\(s\)/.test(L(2)) && /\\mathbf\{I\}/.test(L(2)));
    check('Global basis: ilk adım u_h = Σ a_j N_j, 2. adım N_j tanımı, 3. adım süreklilik (kâğıttan)',
      /\\sum_\{j=1\}/.test(join(at('global', 0))) && /\\begin\{cases\}/.test(join(at('global', 1))) && /continuity/.test(join(at('global', 2))));
    // çizimdeki referans eleman kâğıdın gösterimiyle aynı: ξ yok, s var; eksen değerleri s_k = (k−1)h/p
    const el = at('element', 1, { element: 2 });
    check('Element çizimi ve notu kâğıdın s gösterimini kullanır (ξ yok)', />s<\/text>/.test(el.svg) && !/ξ/.test(el.svg) && /0 ≤ s ≤ h/.test(el.note) && !/ξ/.test(el.note));
    const noSheet = (id, k) => SIMA.renderFrame(b.data, sub, sub.find((q) => q.id === id).offset + k, { width: 700, sheet: null });
    check('kâğıt yoksa yerleşik formül de s gösterimindedir (ξ yok) ve kaynak null',
      ['domain', 'mesh', 'element', 'local', 'global'].every((id) => { const f = noSheet(id, 1); return f.formulaSource === null && f.formulas.length && !/xi/.test(f.formulas.join('')); })
      && !/ξ/.test(sub.find((q) => q.id === 'element').text));
    // son ağdaki s düğüm değerleri çekirdekteki ξ'den türetilmiş s_k = h(1+ξ_k)/2 ile aynı
    const h = b.data.discretization.h, ref = b.data.basis.referenceNodes;
    check('çizimdeki s düğüm etiketleri (k−1)h/p ile aynı', ref.every((t, kk) => Math.abs(h * (1 + t) / 2 - kk * h / c.p) < 1e-12));
  }
  // sayısal yedek kâğıt: kâğıtta olmayan blok için yedek formüle düşer, çökmez
  {
    const c = CASES[0], secN = parse(sheetFor(c.pr, c.n, c.p, false).text), b = mk(c.pr, c.n, c.p), sub = SIMA.stagesFor(SIMA.buildStages(b.data), 'approx');
    let ok = true;
    for (const st of sub) for (let k = 0; k < st.steps; k++) { const f = SIMA.renderFrame(b.data, sub, st.offset + k, { width: 700, sheet: secN }); if (!f.formulas.length || f.formulas.some((t) => /undefined|NaN/.test(t))) ok = false; }
    check('sayısal yedek kâğıtta panel 1 çökmez; her adımın formülü var (kâğıttan ya da yedek)', ok);
  }
  check('sayfa: kâğıt her iki panelde kullanılır, ikisi de görününce istenir, kaynak satırı iki panelde', /sheet: currentSimSheet\(\),/.test(html) && /io\.observe\(\$\('sim-panel'\)\); io\.observe\(\$\('sol-panel'\)\)/.test(html)
    && !/pk !== 'sol' \? ''/.test(html));
}

// 15. Açık / koyu görünüm: sağ üstte Auto-Light-Dark seçici, tercih kalıcı, ilk boyamadan önce uygulanır, durum tek nesnede
{
  const head = html.slice(0, html.indexOf('<body>'));
  check('görünüm seçici başlığın sağında, üç seçenek (auto, light, dark), erişilebilir radiogroup',
    /<header class="top">\s*<div class="top-row">[\s\S]*?<div class="theme" role="radiogroup" aria-label="Appearance">/.test(html) && ['auto', 'light', 'dark'].every((k) => new RegExp(`data-theme-choice="${k}"[^>]*aria-label=`).test(html)));
  check('tercih ilk boyamadan önce (head içinde) uygulanır; yalnız light/dark geçerli, kayıt hatası yutulur', /localStorage\.getItem\('fem-theme'\)[\s\S]*?setAttribute\('data-theme', t\)/.test(head) && /try \{ var t/.test(head));
  check('durum tek nesnede: state.theme, render() applyTheme çağırır; auto özniteliği kaldırır; sistem değişimi izlenir',
    /theme: \(\(\) => \{ try \{ const t = localStorage\.getItem/.test(html) && /function render\(\) \{\s*applyTheme\(\);/.test(html) && /removeAttribute\('data-theme'\)/.test(html) && /addEventListener\('change', \(\) => \{ if \(state\.theme === 'auto'\) render\(\); \}\)/.test(html));
  check('açık ve koyu belirteç kümeleri aynı anahtarlara sahip (biri eksik kalmasın)', (() => {
    const keys = (re) => [...(html.match(re) || [''])[0].matchAll(/(--[a-z0-9-]+):/g)].map((m) => m[1]).sort().join();
    const light = keys(/:root \{[\s\S]*?\n\}/), dark = keys(/:root\[data-theme="dark"\] \{[\s\S]*?\n\}/), mq = keys(/@media \(prefers-color-scheme: dark\) \{[\s\S]*?\n\}\n:root\[data-theme/);
    const L = new Set(light.split(',')), D = new Set(dark.split(','));
    return dark === mq && [...D].every((k) => L.has(k));
  })());
  check('yazdırma koyu seçiliyken de beyaz zemin + siyah yazı', /@media print \{[\s\S]*?:root, :root\[data-theme\] \{ --ink: #000;[\s\S]*?color-scheme: light/.test(html));
}

// 16. Paylaşım öncesi hata düzeltmeleri geri gelmesin: İngilizce kâğıtta Türkçe yok, uç etiketleri sığar, matris başlığı ayrı satırda
{
  const r = spawnSync('python3', ['fem_sheet.py'], { cwd: __dirname, input: JSON.stringify({ ...PROBLEMS[0].inputs, n: 3, degree: 2, symbolic: true }), encoding: 'utf8', timeout: 240000 });
  const sheet = JSON.parse(r.stdout).text, tr = sheet.match(/[çğışÇĞİŞ]|\b(mertebe|interpolasyon|matrisi|sonlu|uyumluluk)\b/g);
  check('İngilizce çözüm kâğıdında Türkçe sözcük yok', sheet.includes('interpolation matrix') && !tr, String(tr && tr.slice(0, 5)));
  check('yakınsama grafiğinin sağ boşluğu en uzun uç etiketine göre (sabit 62 px "theory h^3" yazısını kesiyordu)', /const endChars = Math\.max\(0, \.\.\.cfg\.series\.map/.test(html) && /Math\.max\(62, Math\.ceil\(endChars \* 6\.7\) \+ 14\)/.test(html));
  const b = (() => { const pr = PROBLEMS[0], problem = F.buildProblem(pr.inputs), mesh = F.makeMesh(pr.inputs.x0, pr.inputs.x1, 4, 2), res = F.METHODS.galerkin.solve(problem, mesh);
    return SIMA.buildSimulationData({ problem, mesh, c: res.c, xs: F.sampleXs(mesh), solved: { K: res.K, F: res.F, condition: res.diagnostics.condition }, reference: null, tex: { equation: 'a', conditions: 'b' } }); })();
  const stages = SIMA.stagesFor(SIMA.buildStages(b), 'solution'), st = stages.find((q) => q.id === 'elementK');
  const svg = SIMA.renderFrame(b, stages, st.offset + 2, { width: 700 }).svg, yOf = (re) => Number((svg.match(re) || [])[1]);
  const title = yOf(/<text[^>]*\sy="([\d.]+)"[^>]*>K⁽\d+⁾<\/text>/), column = Math.min(...[...svg.matchAll(/<text[^>]*\sy="([\d.]+)"[^>]*class="sim-node"[^>]*>2<\/text>/g)].map((m) => Number(m[1])));   // "2" hem satır hem sütun numarasıdır; sütun olanı en üsttekidir
  check('eleman matrisi başlığı sütun numaralarının üstünde ayrı satırda (≥ 13 px)', column - title >= 13, `başlık y=${title}, sütun y=${column}`);
}

console.log(`\n${failures ? failures + ' BAŞARISIZ' : 'Tümü geçti'}`);
process.exit(failures ? 1 : 0);
