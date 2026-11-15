// Tarayıcı çekirdeğinin (index.html içindeki CORE bloğu) testleri:  node tests.js
// Aynı dosya, Python testleri için sayısal sonuçları da dışarı verir:  node tests.js --solve '<JSON>'
const fs = require('fs');
const html = fs.readFileSync(__dirname + '/index.html', 'utf8');
const core = html.split('/*CORE-START*/')[1].split('/*CORE-END*/')[0];
const api = new Function(core + '; return { parseExpression, gaussLegendre, shapeFunctions, shapeDerivatives, makeMesh, buildProblem, METHODS, evaluate, referenceSolution, errorNorms, maxErrorRefined, evaluateDerivative, convergenceStudy, observedRate, equationTex, conditionsTex, QUADRATURE_WARNING, leadingDegeneracy, parseExact, analyzeSystem, conditionEstimate, solveLinear, Diagnosis, RANK_TOL, COMPAT_TOL };')();

const D = (value) => ({ type: 'dirichlet', value: String(value) }), Nm = (value) => ({ type: 'neumann', value: String(value) }), NONE = { type: 'none', value: '' };
const problem = (alpha, beta, gamma, f, left, right, x0 = 0, x1 = 1) => api.buildProblem({ alpha, beta, gamma, f, x0, x1, left, right });
const derivativeAt = (s, x) => api.evaluateDerivative(s.mesh, s.c, x);
const solve = (pr, n, p) => { const mesh = api.makeMesh(pr.x0, pr.x1, n, p); return { mesh, ...api.METHODS.galerkin.solve(pr, mesh) }; };

if (process.argv[2] === '--solve-diag') {      // Python testleri için: 1. mertebe uyumluluk teşhisi
  const input = JSON.parse(process.argv[3]);
  console.log(JSON.stringify(solve(api.buildProblem(input), input.n, input.degree).diagnostics.firstOrder));
  process.exit(0);
}
if (process.argv[2] === '--solve') {           // Python tarafı için: verilen problemi çöz, düğüm değerlerini yaz
  const input = JSON.parse(process.argv[3]);
  console.log(JSON.stringify(solve(api.buildProblem(input), input.n, input.degree).c));
  process.exit(0);
}

let fail = 0;
const check = (name, ok, info = '') => { console.log((ok ? 'OK   ' : 'FAIL ') + name, info); if (!ok) fail++; };
const throws = (fn) => { try { fn(); return false; } catch { return true; } };

// 1. Ayrıştırıcı
for (const [expr, x, want] of [['2x^2 - 3x + 1', 2, 3], ['-x^2', 3, -9], ['sin(pi*x)^2', 0.5, 1], ['1/(1+25*x^2)', 1, 1 / 26], ['e^-x', 1, Math.exp(-1)], ['(x+1)(x-1)', 3, 8], ['1/2x', 4, 2]])
  check('parse ' + expr, Math.abs(api.parseExpression(expr).f(x) - want) < 1e-12);
for (const bad of ['sin(', '2+*3', 'foo(x)', 'x)']) check('rejects ' + bad, throws(() => api.parseExpression(bad)));

// 2. Kuadratür ve şekil fonksiyonları
const g = api.gaussLegendre(12);
check('gauss x^22 exact', Math.abs(g.points.reduce((s, x, i) => s + g.weights[i] * x ** 22, 0) - 2 / 23) < 1e-13);
for (const p of [1, 2, 3, 4]) {
  const N = api.shapeFunctions(p, 0.37), dN = api.shapeDerivatives(p, 0.37), eps = 1e-6;
  const numeric = api.shapeFunctions(p, 0.37 + eps).map((v, k) => (v - api.shapeFunctions(p, 0.37 - eps)[k]) / (2 * eps));
  check(`shape p=${p}: ΣN=1, ΣN'=0, N' doğru`, Math.abs(N.reduce((a, b) => a + b) - 1) < 1e-13 && Math.abs(dN.reduce((a, b) => a + b)) < 1e-12 && dN.every((v, k) => Math.abs(v - numeric[k]) < 1e-7));
}

// 3. Elle doğrulanabilir durum: u'' + u = x, u(0)=u(1)=0, n=2, p=1 → a₂ = −3/44 (kâğıtta kesirlerle çıkan sonuç)
check("hand case a2 = -3/44", Math.abs(solve(problem('1', '0', '1', 'x', D(0), D(0)), 2, 1).c[1] + 3 / 44) < 1e-14);

// 4. Tam çözümü bilinen problemler: hata, n ile kuramsal hızda küçülmeli
const cases = [
  ["u''+u=x, D/D", problem('1', '0', '1', 'x', D(0), D(0)), (x) => x - Math.sin(x) / Math.sin(1)],
  ["-u''=sin(pi x), D/D", problem('-1', '0', '0', 'sin(pi*x)', D(0), D(0)), (x) => Math.sin(Math.PI * x) / Math.PI ** 2],
  ["u''-4u=0, D/N", problem('1', '0', '-4', '0', D(0), Nm(1)), (x) => Math.sinh(2 * x) / (2 * Math.cosh(2))],
  ["((1+x)u')'=x, D/N", problem('1+x', '0', '0', 'x', D(0), Nm(0)), (x) => x * x / 4 - x / 2],
  ["-0.05u''+u'=1, D/D", problem('-0.05', '1', '0', '1', D(0), D(0)), (x) => x - (Math.exp(20 * x) - 1) / (Math.exp(20) - 1)],
  ["u''=2, N/D (u=x^2-2x)", problem('1', '0', '0', '2', Nm(-2), D(-1)), (x) => x * x - 2 * x]
];
for (const [name, pr, exact] of cases) {
  for (const p of [1, 2, 3, 4]) {
    const errors = [8, 16].map((n) => { const s = solve(pr, n, p); return api.errorNorms(exact, s.mesh, s.c).l2; });
    const rate = errors[1] < 1e-12 ? null : Math.log2(errors[0] / errors[1]);
    check(`${name} p=${p}`, rate === null || rate > p + 0.6, rate === null ? 'tam' : `hız ${rate.toFixed(2)}, hata ${errors[1].toExponential(1)}`);
  }
  const ref = api.referenceSolution(pr, 'galerkin');
  const gap = Math.max(...Array.from({ length: 101 }, (_, i) => Math.abs(ref(i / 100) - exact(i / 100))));
  check(`${name}: referans çözüm tam çözüme yakın`, gap < 1e-8, gap.toExponential(1));
}
// derecesi ≤ p olan tam çözüm aynen bulunmalı
{ const s = solve(cases[5][1], 3, 2); check('polinom çözüm tam bulunur', api.errorNorms(cases[5][2], s.mesh, s.c).max < 1e-12); }

// 5. Birinci mertebe: u' + 2u = 0, u(0)=1 → e^{-2x}
{ const pr = problem('0', '1', '2', '0', D(1), NONE), exact = (x) => Math.exp(-2 * x);
  for (const p of [1, 2, 3]) {
    const errors = [8, 16].map((n) => { const s = solve(pr, n, p); return api.errorNorms(exact, s.mesh, s.c).l2; });
    check(`1. mertebe p=${p}: hata küçülüyor`, errors[1] < errors[0] / 3 && errors[1] < 1e-2, errors.map((v) => v.toExponential(1)).join(' → '));
  } }

// 5b. İkinci mertebe başlangıç/son değer problemi (aynı uçta u ve u'): L² hızı p + 1, analitik çözüme karşı
const B = (u, s) => ({ type: 'both', value: String(u), slope: String(s) });
const ivpCases = [
  ["IVP u''+u=0, u(0)=0, u'(0)=1", problem('1', '0', '1', '0', B(0, 1), NONE), Math.sin],
  ["IVP u''+u=0 uzun aralık [0,10]", problem('1', '0', '1', '0', B(0, 1), NONE, 0, 10), Math.sin],
  ["IVP sönümlü u''+0.4u'+4u=0", problem('1', '0.4', '4', '0', B(1, 0), NONE, 0, 5),
    (x) => { const w = 3 * Math.sqrt(11) / 5; return Math.exp(-x / 5) * (Math.cos(w * x) + Math.sqrt(11) / 33 * Math.sin(w * x)); }],
  ["IVP ((1+x)u')'=0", problem('1+x', '0', '0', '0', B(0, 1), NONE), (x) => Math.log(1 + x)],
  ["son değer u''-u=0, u(1)=1, u'(1)=0", problem('1', '0', '-1', '0', NONE, B(1, 0)), (x) => Math.cosh(x - 1)]
];
for (const [name, pr, exact] of ivpCases) {
  for (const p of [1, 2, 3, 4]) {
    const errors = [8, 16].map((n) => { const s = solve(pr, n, p); return api.errorNorms(exact, s.mesh, s.c).l2; });
    const rate = errors[1] < 1e-11 ? null : Math.log2(errors[0] / errors[1]);
    check(`${name} p=${p}`, rate === null || rate > p + 0.6, rate === null ? `yuvarlama tabanı, hata ${errors[1].toExponential(1)}` : `hız ${rate.toFixed(2)}, hata ${errors[1].toExponential(1)}`);
  }
}
{ const s = solve(problem('1', '0', '0', '-9.81', B(0, 20), NONE, 0, 4), 3, 2);   // atış: u = 20x − 4.905x² (p ≥ 2'de tam)
  check('IVP atış problemi tam', api.errorNorms((x) => 20 * x - 4.905 * x * x, s.mesh, s.c).max < 1e-10); }
check('IVP: öbür uca koşul reddedilir', throws(() => problem('1', '0', '1', '0', B(0, 1), D(0))));
check('IVP: iki uçta u ve u\' reddedilir', throws(() => problem('1', '0', '1', '0', B(0, 1), B(0, 1))));
check('IVP: 1. mertebede reddedilir', throws(() => problem('0', '1', '1', '0', B(0, 1), NONE)));
check('IVP: α(x₀)=0 reddedilir', throws(() => problem('x', '0', '1', '0', B(0, 1), NONE)));

// 5c. Kıyas problemi u' = 2x, u(0) = 0 (tam çözüm x²), n = 2
{ const pr = problem('0', '1', '0', '2x', D(0), NONE), exact = (x) => x * x, dexact = (x) => 2 * x;
  const s1 = solve(pr, 2, 1);
  check("kıyas P1: a = [0, 1/6, 1]", [0, 1 / 6, 1].every((v, i) => Math.abs(s1.c[i] - v) < 1e-14), s1.c.join(', '));
  check("kıyas P1: K^(1) ders kitabı işaretiyle [[-1/2, 1/2], [-1/2, 1/2]] (1. mertebe, kısmi integrasyon yok)",
    Math.abs(s1.K[0][0] + 0.5) < 1e-14 && Math.abs(s1.K[0][1] - 0.5) < 1e-14 && Math.abs(s1.K[1][0] + 0.5) < 1e-14 && Math.abs(s1.F[0] - 1 / 12) < 1e-14);
  const n1 = api.errorNorms(exact, s1.mesh, s1.c, dexact);
  check('kıyas P1: düğüm hatası 1/12, L² = 0.0304290, H¹ = 1/3', Math.abs(n1.nodal - 1 / 12) < 1e-14 && Math.abs(n1.l2 - Math.sqrt(0.000925925925925926)) < 1e-12 && Math.abs(n1.h1 - 1 / 3) < 1e-12,
    `nodal ${n1.nodal}, l2 ${n1.l2}, h1 ${n1.h1}`);
  check('kıyas P1: u_h = x/3 ve 5x/3 − 2/3', [0.1, 0.3, 0.6, 0.9].every((x) => Math.abs(api.evaluate(s1.mesh, s1.c, x) - (x <= 0.5 ? x / 3 : 5 * x / 3 - 2 / 3)) < 1e-14));
  for (const p of [2, 3]) {
    const sp = solve(pr, 2, p), nm = api.errorNorms(exact, sp.mesh, sp.c, dexact);
    check(`kıyas P${p}: N = ${2 * p + 1}, x² makine hassasiyetinde`, sp.c.length === 2 * p + 1 && nm.l2 < 1e-14 && nm.h1 < 1e-13, `L² ${nm.l2.toExponential(1)}, H¹ ${nm.h1.toExponential(1)}`);
  } }

// 5d. Şekil fonksiyonları: N_i(x_j) = δ_ij (eşit aralıklı düğümler: P2 orta nokta, P3 üçte birler)
for (const p of [1, 2, 3, 4]) {
  const nodes = Array.from({ length: p + 1 }, (_, k) => -1 + 2 * k / p);
  check(`shape p=${p}: N_i(x_j) = δ_ij`, nodes.every((xj, j) => api.shapeFunctions(p, xj).every((v, i) => Math.abs(v - (i === j ? 1 : 0)) < 1e-14)));
}

// 5e. H¹ yarı-norm hızı p (2. mertebe sınır değer problemi) ve gözlenen hız gerçek h ile
{ const pr = cases[1][1], exact = cases[1][2], dexact = (x) => Math.cos(Math.PI * x) / Math.PI;
  for (const p of [1, 2, 3]) {
    const e = [8, 16].map((n) => { const s = solve(pr, n, p); return api.errorNorms(exact, s.mesh, s.c, dexact).h1; });
    const rate = Math.log2(e[0] / e[1]);
    check(`H¹ hızı p=${p} ≈ ${p}`, Math.abs(rate - p) < 0.15, `hız ${rate.toFixed(2)}`);
  } }
check('observedRate: gerçek h ile (h oranı 3)', Math.abs(api.observedRate([{ h: 0.3, E: 9e-2 }, { h: 0.1, E: 1e-2 }]).rate - 2) < 1e-12);
check('observedRate: tabandaki noktalar kullanılmaz, hız tanımsız (0 değil)', api.observedRate([{ h: 0.5, E: 1e-15 }, { h: 0.25, E: 2e-16 }], () => 1e-12) === null);
check('observedRate: taban üstündeki son çift seçilir', Math.abs(api.observedRate([{ h: 1, E: 1e-2 }, { h: 0.5, E: 2.5e-3 }, { h: 0.25, E: 1e-15 }], () => 1e-12).rate - 2) < 1e-12);
check("α(x₀)=0 ucunda u' koşulu reddedilir", throws(() => problem('x', '0', '0', '0', Nm(1), D(1))));

// 5f. E_max: örnek noktaları arasında kalan tepe bulunur (u = x e^{-3x}, u_h = 0 → max x = 1/3'te e^{-1}/3; 1/3 örnek noktası değil)
{ const mesh = api.makeMesh(0, 1, 1, 1), m = api.maxErrorRefined((x) => x * Math.exp(-3 * x), mesh, [0, 0]);
  check('E_max yerel arama: kesin tepe 1e-12 içinde', Math.abs(m.value - Math.exp(-1) / 3) < 1e-12 && Math.abs(m.x - 1 / 3) < 1e-5, `${m.value} @ ${m.x}`); }

// 5g. IVP kararlılık sınırı (p = 1, u'' + ω²u = 0): üç terimli yinelemenin kökleri ωh < 2√3 iken birim çemberde.
//     Sınırın altında genlik korunur, üstünde üstel büyür. Bu bir yöntem sınırıdır; test sınırın varlığını belgeler.
{ const amp = (wh) => { const L = 64 * wh, pr = problem('1', '0', '1', '0', B(0, 1), NONE, 0, L), s = solve(pr, 64, 1);
    return Math.max(...s.c.slice(-8).map(Math.abs)); };
  check('IVP p=1: ωh = 3 (< 2√3) genlik sınırlı', amp(3) < 3, amp(3).toFixed(2));
  check('IVP p=1: ωh = 3.6 (> 2√3) üstel büyüme', amp(3.6) > 1e6, amp(3.6).toExponential(1)); }
// Kaba ağdaki IVP tekilliği "çözüm tek değil" diye değil, ağ uyarısıyla bildirilir (p = 2, ω²h² = 10 ayrık rezonansı yakını)
{ let message = '';
  try { solve(problem('1', '0', '1', '0', B(0, 1), NONE, 0, 64 * 8), 64, 1); } catch (err) { message = err.message; }
  check('IVP kaba ağ: doğru mesaj (tekil değil, ağ kaba)', message.includes('ağ bu çözüm için fazla kaba') && !message.includes('En az bir uçta'), message.slice(0, 60)); }
// u'(x₀) koşulu zayıf (doğal) biçimde girer: u_h'(x₀) ≠ s₀ ama ağ inceldikçe yaklaşır
{ const pr = ivpCases[0][1], g = [8, 32].map((n) => { const s = solve(pr, n, 2); return Math.abs(derivativeAt(s, 0) - 1); });
  check("IVP: |u_h'(x₀) − u'(x₀)| ağla küçülür (zayıf koşul)", g[0] > 1e-4 && g[1] < g[0] / 8, g.map((v) => v.toExponential(1)).join(' → ')); }

// 5h. Ayrık sistem analizi (SVD): rank, uyumluluk (geriye hata), ölçekten ve boyuttan bağımsız karar
{ const scaled = (A, s) => A.map((r) => r.map((v) => v * s));
  for (const s of [1e-8, 1, 1e8]) {
    const A = scaled([[1, 1], [1, 1]], s);
    const ok = api.analyzeSystem(A, [s, s]), bad = api.analyzeSystem(A, [s, 2 * s]), I = api.analyzeSystem(scaled([[1, 0], [0, 1]], s), [s, s]);
    check(`analyzeSystem ölçek ${s}: rank 1 uyumlu / rank 1 uyumsuz / rank 2`, ok.rank === 1 && ok.compatible && bad.rank === 1 && !bad.compatible && I.rank === 2,
      `η ${ok.eta.toExponential(1)} / ${bad.eta.toExponential(1)}`);
  }
  // Boyut: 1-B Laplace (Neumann) matrisi, rank N−1, sabit boş uzay; uyumlu ve uyumsuz sağ taraf
  for (const N of [3, 17, 65, 257]) {
    const A = Array.from({ length: N }, (_, i) => Array.from({ length: N }, (_, j) => (i === j ? ((i === 0 || i === N - 1) ? 1 : 2) : (Math.abs(i - j) === 1 ? -1 : 0))));
    const comp = Array.from({ length: N }, (_, i) => Math.sin(i + 1)), mean = comp.reduce((a, b) => a + b) / N;
    const r1 = api.analyzeSystem(A, comp.map((v) => v - mean)), r2 = api.analyzeSystem(A, comp);
    check(`analyzeSystem N=${N}: rank N−1, uyumlu ↔ uyumsuz ayrımı`, r1.rank === N - 1 && r1.compatible && !r2.compatible, `η ${r1.eta.toExponential(1)} / ${r2.eta.toExponential(1)}`);
  }
  check('analyzeSystem: NaN → sonlu değil', api.analyzeSystem([[1, NaN], [0, 1]], [1, 1]).finite === false);
}
// Koşul sayısı tahmini (κ₁): bilinen değerlere karşı (Hilbert 6: κ₁ = 2.907e7; köşegen: 1e6)
{ const H = Array.from({ length: 6 }, (_, i) => Array.from({ length: 6 }, (_, j) => 1 / (i + j + 1)));
  const lu = {}; api.solveLinear(H, new Array(6).fill(1), lu); const k = api.conditionEstimate(H, lu);
  check('conditionEstimate Hilbert(6) 3 kat içinde', k > 2.907e7 / 3 && k <= 2.907e7 * 1.0001, k.toExponential(3));
  const Dg = [[1, 0], [0, 1e-6]], lu2 = {}; api.solveLinear(Dg, [1, 1], lu2);
  check('conditionEstimate köşegen = 1e6', Math.abs(api.conditionEstimate(Dg, lu2) / 1e6 - 1) < 1e-9); }
// Kategoriler: ayrık tekillik nedeniyle ayrıştırılır; genel "tekil" mesajı yok
{ const cat = (fn) => { try { fn(); return 'yok'; } catch (e) { return e.category || 'kategorisiz'; } };
  check("K1: u''=1, iki uç doğal (∫f ≠ 0, uyumsuz)", cat(() => solve(problem('1', '0', '0', '1', NONE, NONE), 4, 1)) === 'K1');
  check("K2: u''=1, u'(0)=0, u'(1)=1 (uyumlu, sabit boş uzay)", cat(() => solve(problem('1', '0', '0', '1', Nm(0), Nm(1)), 4, 2)) === 'K2');
  check("K2: u''=2x−1, N/N (uyumlu)", cat(() => solve(problem('1', '0', '0', '2x-1', Nm(0), Nm(0)), 4, 2)) === 'K2');
  check('K3: IVP kaba ağ (ωh > 2√3)', cat(() => solve(problem('1', '0', '1', '0', B(0, 1), NONE, 0, 64 * 8), 64, 1)) === 'K3');
  check('K1/K2: işaret değiştiren α, n çift (eski genel mesaj yok)', ['K1', 'K2', 'K3'].includes(cat(() => solve(problem('x-0.5', '0', '0', '1', D(0), D(1)), 4, 1))));
  check("kaynakta genel 'Sistem tekil: çözüm tek değil' mesajı kalmadı", !html.includes('Sistem tekil: çözüm tek değil'));
  check('koşul sayısı her çözümde raporlanır', solve(cases[0][1], 8, 2).diagnostics.condition > 1); }

// 5i. Dejenere katsayı tespiti (madde 2a): uyarı düzeyi, ret değil
{ const kind = (a, x0 = 0, x1 = 1) => problem(a, '0', '0', '1', D(0), D(1), x0, x1).degeneracy;
  check('α=1: düzgün', kind('1').kind === 'uniform');
  const e = kind('x'); check('α=x: dejenere, x*=0 uçta', e.kind === 'degenerate' && Math.abs(e.points[0].x) < 1e-9);
  const sc = kind('x-0.5'); check('α=x−0.5: işaret değiştiren, x*=0.5', sc.kind === 'sign_changing' && Math.abs(sc.points[0].x - 0.5) < 1e-9);
  const touch = kind('(x-0.3)^2'); check('α=(x−0.3)²: değme sıfırı (işaret değiştirmez) bulunur', touch.kind === 'degenerate' && Math.abs(touch.points[0].x - 0.3) < 1e-6, String(touch.points[0]?.x));
  const sn = kind('sin(x)', 0, 4); check('α=sin x [0,4]: 0 ve π', sn.kind === 'sign_changing' && sn.points.length === 2 && Math.abs(sn.points[1].x - Math.PI) < 1e-9);
  check('1. mertebe β=x−0.5: işaret değiştiren', problem('0', 'x-0.5', '0', '1', D(0), NONE).degeneracy.kind === 'sign_changing');
  check('dejenerelik ret değil: α=x−0.5, n tek çözülür', !throws(() => solve(problem('x-0.5', '0', '0', '1', D(0), D(1)), 5, 1))); }

// 5j. Analitik çözümün tarayıcıda bağımsız koşul kontrolü (madde 1): sunucudan gelen aday koşulu sağlamıyorsa kullanılmaz
{ const payload = (u, du, f, df) => ({ u, du, tex: '', samples: [0, 0.25, 0.5, 0.75, 1].map((x) => [x, f(x), df(x)]) });
  const bad = api.parseExact(payload('x-1', '1', (x) => x - 1, () => 1), problem('x', '0', '0', '1', D(0), D(0)));
  check('parseExact: u(0)=−1 ≠ 0 adayı reddedilir (KR-1)', bad && bad.rejected && bad.rejected.startsWith('u(0)'), bad && bad.rejected);
  const good = api.parseExact(payload('x-1', '1', (x) => x - 1, () => 1), problem('x', '0', '0', '1', NONE, D(0)));
  check('parseExact: dejenere serbest uçta akı sınanmaz, aday kabul', good && !good.rejected && typeof good.u === 'function');
  const nat = api.parseExact(payload('x', '1', (x) => x, () => 1), problem('1', '0', '0', '0', NONE, D(1)));
  check('parseExact: doğal koşul α u′(0) = 1 ≠ 0 reddedilir', nat && nat.rejected && nat.rejected.includes('doğal'));
  const ivp = api.parseExact(payload('2*sin(x)', '2*cos(x)', (x) => 2 * Math.sin(x), (x) => 2 * Math.cos(x)), problem('1', '0', '1', '0', B(0, 1), NONE));
  check('parseExact: IVP eğimi u′(0) = 2 ≠ 1 reddedilir', ivp && ivp.rejected && ivp.rejected.startsWith('u′(0)')); }

// 5k. Birinci mertebe uyumluluk (madde 3): atılan satır kalıntısı ve sol sıfır uzayı; yᵀb = −r özdeşliği
{ const fo = (f, n, p, g = '0', x1 = 2) => solve(problem('0', '1', g, f, D(g === '0' ? 0 : 1), NONE, 0, x1), n, p).diagnostics.firstOrder;
  const hw = fo('2x', 2, 1);
  check("u'=2x [0,2] n=2 P1 (ödev): r = 0, uyumlu", Math.abs(hw.residual) < 1e-14 && hw.compatible && hw.row === 1, `r ${hw.residual}`);
  const cubic = fo('3x^2', 2, 1);
  check("u'=3x² [0,2] n=2 P1: r = −1, uyumsuz", Math.abs(cubic.residual + 1) < 1e-13 && !cubic.compatible, `r ${cubic.residual}`);
  for (const [f, n, p, g] of [['2x', 2, 1, '0'], ['3x^2', 8, 1, '0'], ['0', 8, 1, '2'], ['3x^2', 4, 3, '0']]) {
    const r = fo(f, n, p, g);
    check(`yᵀb = −r (iki ayrı hesap) f=${f} γ=${g} n=${n} p=${p}`, Math.abs(r.yb + r.residual) <= 1e-12 * Math.max(1, Math.abs(r.residual)), `${r.yb} / ${r.residual}`);
  }
  check('2. mertebede firstOrder yok', solve(cases[0][1], 4, 1).diagnostics.firstOrder === null); }

// 5l. V0.1 son kontroller: işaret değiştiren α ve referans ağ; toleransların ölçek/boyut davranışı (FEM matrisleriyle)
{ const sign = problem('x-0.5', '0', '0', '1', D(0), D(1));
  check('işaret değiştiren α, kullanıcı ağı n=5 (x=0.5 düğüm değil): çözülür', !throws(() => solve(sign, 5, 1)));
  let ref = null; try { api.referenceSolution(sign, 'galerkin'); } catch (e) { ref = e; }
  check('aynı problemde referans ağ n=128 (x=0.5 düğüm): kategorili teşhis, sayı üretilmez', ref && ['K1', 'K2', 'K3'].includes(ref.category), ref && ref.message.slice(0, 60));
  let own = null; try { solve(sign, 4, 1); } catch (e) { own = e; }
  check('kullanıcı ağı n=4 (x=0.5 düğüm): kategorili teşhis', own && ['K1', 'K2', 'K3'].includes(own.category), own && own.category);
  // Ölçek: α ve f birlikte 1e-8…1e8 ölçeklenir; saf Neumann uyumlu → K2, uyumsuz → K1 (karar ölçekten bağımsız olmalı)
  for (const k of [-8, 0, 8]) {
    const sc = `10^${k}`, cat = (f) => { try { solve(problem(sc, '0', '0', f, Nm(0), Nm(0)), 8, 2); return 'yok'; } catch (e) { return e.category; } };
    check(`ölçek ${sc}: uyumlu K2, uyumsuz K1`, cat(`${sc}*(2x-1)`) === 'K2' && cat(sc) === 'K1');
  }
  check('bilimsel gösterim 1e-8 reddedilir (sessizce 1·e − 8 okunmaz)', throws(() => api.parseExpression('1e-8')) && throws(() => api.parseExpression('2.5E3')) && !throws(() => api.parseExpression('2exp(x)')) && !throws(() => api.parseExpression('3e')));
  // Boyut: N = n·p + 1 = 3 … 257
  const sizes = [[2, 1], [8, 2], [32, 4], [64, 4]].map(([n, p]) => { try { solve(problem('1', '0', '0', '2x-1', Nm(0), Nm(0)), n, p); return 'yok'; } catch (e) { return e.category; } });
  check('boyut N=3…257: uyumlu saf Neumann her boyutta K2', sizes.every((c) => c === 'K2'), sizes.join(','));
  // κ tahmini ölçekten bağımsız (göreli büyüklük)
  const k1 = solve(problem('1', '0', '1', 'x', D(0), D(0)), 16, 3).diagnostics.condition, k2 = solve(problem('10^8', '0', '10^8', '10^8*x', D(0), D(0)), 16, 3).diagnostics.condition;
  check('κ̂ ölçekten bağımsız (α, γ, f × 1e8)', Math.abs(k1 / k2 - 1) < 1e-6, `${k1.toExponential(3)} / ${k2.toExponential(3)}`); }

// 6. Geçersiz ya da çözümü tek olmayan girdiler açık hata vermeli
check('tekil: u\'\'=1, koşul yok', throws(() => solve(problem('1', '0', '0', '1', NONE, NONE), 4, 1)));
check('tekil: u\'\'=1, iki uç da u\'', throws(() => solve(problem('1', '0', '0', '1', Nm(0), Nm(1)), 4, 2)));
check('koşulsuz ama γ≠0 çözülür', !throws(() => solve(problem('1', '0', '-1', '1', NONE, NONE), 4, 1)));
check('1. mertebe koşulsuz reddedilir', throws(() => problem('0', '1', '1', '0', NONE, NONE)));
check('1. mertebe u\' koşulu reddedilir', throws(() => problem('0', '1', '1', '0', Nm(0), NONE)));
check('α=β=0 reddedilir', throws(() => problem('0', '0', '1', 'x', D(0), D(0))));
check('aralıkta tanımsız katsayı reddedilir', throws(() => problem('1', '0', '1/x', '0', D(0), D(0))));
check('x₀ ≥ x₁ reddedilir', throws(() => problem('1', '0', '0', '1', D(0), D(0), 1, 0)));

// 7. Kuadratür uyarısı: düzgün katsayıda yok, köşeli yükte var
check('kuadratür: düzgün ≈ 0', solve(cases[1][1], 4, 2).quadratureError < 1e-12);
check('kuadratür: abs(x-0.3) işaretlenir', solve(problem('1', '0', '0', 'abs(x-0.3)', D(0), D(0)), 3, 2).quadratureError > api.QUADRATURE_WARNING);

// 8. Önizleme metni
check('equationTex', api.equationTex(cases[0][1]) === "u'' + u = x" && api.equationTex(cases[4][1]).startsWith("-0.05\\,u'' + u'"), api.equationTex(cases[3][1]));
process.exit(fail ? 1 : 0);
