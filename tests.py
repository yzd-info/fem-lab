#!/usr/bin/env python3
"""Çözüm kâğıdı testleri:  python3 tests.py

Üç bağımsız kaynak birbirine karşı sınanır:
  1. fem_sheet.py   — SymPy ile sembolik (kâğıttaki sayılar)
  2. index.html     — JavaScript ile sayısal (grafiklerdeki sayılar; `node tests.js --solve` üzerinden)
  3. elle bilinen sonuçlar ve analitik çözümler
pdflatex kuruluysa her kâğıdın .tex çıktısı derlenir.
"""
import json
import shutil
import subprocess
import sys
import tempfile
import time
from pathlib import Path

import sympy as sp

import fem_sheet
import fem_verify

HERE = Path(__file__).resolve().parent
D = lambda v: {'type': 'dirichlet', 'value': str(v)}
Nm = lambda v: {'type': 'neumann', 'value': str(v)}
NONE = {'type': 'none', 'value': ''}
B = lambda u, s: {'type': 'both', 'value': str(u), 'slope': str(s)}   # aynı uçta u ve u' (başlangıç değer problemi)
failures = 0


def check(name, ok, info=''):
    global failures
    print(('OK   ' if ok else 'FAIL ') + name, info)
    failures += not ok


def problem(alpha, beta, gamma, f, left, right, n, degree, x0=0, x1=1):
    return dict(alpha=alpha, beta=beta, gamma=gamma, f=f, x0=x0, x1=x1, left=left, right=right, n=n, degree=degree)


def javascript_solution(data):
    out = subprocess.run(['node', str(HERE / 'tests.js'), '--solve', json.dumps(data)], capture_output=True, text=True, check=True)
    return json.loads(out.stdout)


# 1. Elle bilinen sonuç: u'' + u = x, u(0)=u(1)=0, n=2, p=1  →  K = [[11/6, −25/12, 0], …],  a₂ = −3/44
r = fem_sheet.solve_problem(problem('1', '0', '1', 'x', D(0), D(0), 2, 1))
check('elle: K₁₁ = 11/6, K₁₂ = −25/12', r['K'][0, 0] == sp.Rational(11, 6) and r['K'][0, 1] == sp.Rational(-25, 12))
check('elle: a₂ = −3/44', r['a'][1] == sp.Rational(-3, 44), str(r['a']))

x = fem_sheet.x
# 1b. Kıyas problemi u' = 2x, u(0) = 0 (tam çözüm x²), n = 2: birinci mertebe, ders kitabı formülasyonu
B1 = lambda p: problem('', '1', '', '2x', D(0), NONE, 2, p)
r = fem_sheet.solve_problem(B1(1))
check('kıyas P1: K^(1) = [[-1/2, 1/2], [-1/2, 1/2]], F^(1) = [1/12, 1/6]',
      r['elements'][0]['K'] == sp.Matrix([[-sp.Rational(1, 2), sp.Rational(1, 2)], [-sp.Rational(1, 2), sp.Rational(1, 2)]])
      and list(r['elements'][0]['F']) == [sp.Rational(1, 12), sp.Rational(1, 6)], str(r['elements'][0]['K'].tolist()))
check('kıyas P1: a = [0, 1/6, 1]', r['a'] == [0, sp.Rational(1, 6), 1], str(r['a']))
pieces = [sp.expand(sum(r['a'][e + k] * r['N'][k] for k in range(2)).subs(fem_sheet.s, x - e * r['h'])) for e in range(2)]
check('kıyas P1: u_h = x/3, 5x/3 − 2/3', pieces == [x / 3, 5 * x / 3 - sp.Rational(2, 3)], str(pieces))
sheet = fem_sheet.make_sheet(B1(1))
check('kıyas P1 kâğıt: düğüm hatası 1/12 kesir, L² satırı, 1. mertebe formülasyonu (α ve F_b yok)',
      '\\frac{1}{12}' in sheet['text'] and '\\|u - u_h\\|_{L^2}' in sheet['text'] and '0.030429' in sheet['text']
      and 'kısmi integrasyon gerekmez' in sheet['text'] and 'F_{b,i}' not in sheet['text'] and "N_i'\\,\\alpha" not in sheet['text'])
for p in (2, 3):
    rp = fem_sheet.solve_problem(B1(p))
    nodes = [sp.Rational(i, 2 * p) for i in range(2 * p + 1)]
    check(f'kıyas P{p}: N = {2 * p + 1}, a_i = x_i² tam', len(rp['a']) == 2 * p + 1 and all(a == t ** 2 for a, t in zip(rp['a'], nodes)))
# Şekil fonksiyonları: N_i(x_j) = δ_ij, ΣN = 1, ΣN' = 0
h_ = sp.Rational(1, 3)
for p in (1, 2, 3, 4):
    Nf, dNf = fem_sheet.lagrange_basis(p, h_)
    ok = all(sp.simplify(Nf[i].subs(fem_sheet.s, sp.Rational(j, p) * h_) - (1 if i == j else 0)) == 0 for i in range(p + 1) for j in range(p + 1))
    check(f'şekil p={p}: δ_ij, ΣN = 1, ΣN\' = 0', ok and sp.simplify(sum(Nf) - 1) == 0 and sp.simplify(sum(dNf)) == 0)
# Tarayıcıya giden analitik çözüm metni ayrıştırılabilir biçimde ('**' yok, ondalık yok)
ej = fem_sheet.exact_only(problem('-0.05', '1', '0', '1', D(0), D(0), 1, 1))['exact_js']
check('exact_js: ^ biçiminde, örnek değerler var', ej is not None and '**' not in ej['u'] and len(ej['samples']) == 9, ej and ej['u'])

# 1c. IVP: u'(x₀) zayıf biçimde sağlanır ve kâğıt bunu gösterir; kaba ağ tekilliği doğru mesajla bildirilir
t = fem_sheet.make_sheet(problem('1', '0', '1', '0', B(0, 1), NONE, 2, 1))['text']
check("IVP kâğıt: u_h'(0) = 24/25 ile u'(0) = 1 yan yana (zayıf koşul)", "u_h'(0) = \\frac{24}{25}" in t and "u'(0) = 1" in t)
try:
    fem_sheet.solve_problem(problem('1', '0', '1', '0', B(0, 1), NONE, 64, 1, 0, 512), symbolic=False)
    check('IVP kaba ağ: doğru mesaj', False)
except fem_sheet.UserError as error:
    check('IVP kaba ağ: doğru mesaj', 'fazla kaba' in str(error) and 'En az bir uçta' not in str(error), str(error)[:60])

# 1d. Ayrık sistem teşhis kategorileri (madde 4): tam yolda kesin rank, sayısal yolda SVD + geriye hata
def category(data, symbolic=True):
    try:
        fem_sheet.solve_problem(data, symbolic=symbolic); return 'çözüldü'
    except fem_sheet.UserError as error:
        return getattr(error, 'category', 'kategorisiz')
check("K1 tam yol: u''=1 koşulsuz (∫f ≠ 0)", category(problem('1', '', '', '1', NONE, NONE, 2, 1)) == 'K1')
check("K2 tam yol: u''=1, u'(0)=0, u'(1)=1 (uyumlu)", category(problem('1', '', '', '1', Nm(0), Nm(1), 2, 2)) == 'K2')
check("K2 sayısal yol: u''=2x−1, N/N, n=4 p=2", category(problem('1', '', '', '2x-1', Nm(0), Nm(0), 4, 2)) == 'K2')
check("K1 sayısal yol: u''=1 koşulsuz, n=4 p=2", category(problem('1', '', '', '1', NONE, NONE, 4, 2)) == 'K1')
check('K3: IVP kaba ağ', category(problem('1', '', '1', '0', B(0, 1), NONE, 64, 1, 0, 512), symbolic=False) == 'K3')
check("kaynakta genel 'Sistem tekil: çözüm tek değil' mesajı kalmadı", 'Sistem tekil: çözüm tek değil' not in (HERE / 'fem_sheet.py').read_text(encoding='utf-8'))
try:
    fem_sheet.solve_problem(problem('x-0.5', '', '', '1', D(0), D(1), 4, 1))
except fem_sheet.UserError as error:
    check('K2 önerisi u verilmişken yazılmaz', 'Bir uçta u değeri vermek' not in str(error), str(error)[:70])

# 1e. Sürekli sınıflandırma (madde 2b): yalnız sembolik olarak kanıtlanan K1/K2 kâğıdı engeller
def continuous(data):
    try:
        r = fem_sheet.solve_problem({**data, 'n': 2, 'degree': 1}, classify=True)
        return r['classification']['kind'], r['classification'].get('method'), r['degeneracy']['kind']
    except fem_sheet.Diagnosis as error:
        return error.category, error.level, None
cont = lambda a, b, g, f, l, r, x0=0, x1=1: continuous(dict(alpha=a, beta=b, gamma=g, f=f, x0=x0, x1=x1, left=l, right=r))
check('E1 rezonans, ∫f sin πx ≠ 0 → K1 sürekli', cont('1', '', 'pi^2', 'sin(pi x)', D(0), D(0))[:2] == ('K1', 'sürekli problem'))
check('E1b rezonans, f=0 → K2 sürekli', cont('1', '', 'pi^2', '0', D(0), D(0))[:2] == ('K2', 'sürekli problem'))
check('E1c uyumlu rezonans (sin 2πx) → K2 sürekli', cont('1', '', 'pi^2', 'sin(2 pi x)', D(0), D(0))[:2] == ('K2', 'sürekli problem'))
check('E1d rezonansa yakın → tek, sembolik (büyüklükten çözümsüzlük iddiası yok)', cont('1', '', 'pi^2+1/1000000', 'sin(pi x)', D(0), D(0))[:2] == ('unique', 'sembolik'))
check('E2 α=x, iki uçta u → K1 (sınırlı çözüm varsayımıyla)', cont('x', '', '', '1', D(0), D(0))[:2] == ('K1', 'sürekli problem'))
check('E4 saf Neumann uyumlu → K2; E4b uyumsuz → K1', cont('1', '', '', '2x-1', Nm(0), Nm(0))[0] == 'K2' and cont('1', '', '', '1', Nm(0), Nm(0))[0] == 'K1')
check('T7 değişken katsayı → belirlenemedi (ret yok)', cont('1+x^2', 'sin(x)', 'exp(x)', '1', D(0), D(1))[0] == 'undetermined')
check('işaret değiştiren α → belirlenemedi + uyarı (ret yok)', cont('x-0.5', '', '', '1', D(0), D(1), 0, 1)[0] in ('undetermined', 'K2'))
check('T12 dejenere uç serbest → tek', cont('x', '', '', '1', NONE, D(0))[:1] == ('unique',))
check('1. mertebe IVP → tek (m = 1)', cont('', '1', '', '2x', D(0), NONE, 0, 2)[:2] == ('unique', 'sembolik'))
check('2. mertebe IVP → tek', cont('1', '', '1', '0', B(0, 1), NONE)[:2] == ('unique', 'sembolik'))

# 1f. Analitik çözüm doğrulaması (madde 1): yapay adaylar doğru nedenle reddedilir; benchmark analitikleri doğrulanır
def verdict(a, b, g, f, l, r, cand, x0=0, x1=1):
    res = fem_sheet.solve_problem(problem(a, b, g, f, l, r, 2, 1, x0, x1))
    return fem_verify.verify_exact(res, cand)
X = fem_sheet.x
for name, args, cand, start in [("u''=2, D/D, aday x² (u(1) ihlali)", ('1', '', '', '2', D(0), D(0)), X**2, 'u(1)'),
                                ("u''=2, aday x³ (denklem ihlali)", ('1', '', '', '2', D(0), D(1)), X**3, 'R[u]'),
                                ("u''=0, u'(1)=1, aday 0", ('1', '', '', '0', D(0), Nm(1)), sp.Integer(0), "u′(1)"),
                                ("u''=0, sol verilmedi, aday x (doğal koşul ihlali)", ('1', '', '', '0', NONE, D(1)), X, 'α u′(0)'),
                                ("IVP, aday 2 sin x (eğim ihlali)", ('1', '', '1', '0', B(0, 1), NONE), 2 * sp.sin(X), "u′(0)"),
                                ("1. mertebe, aday x²+1 (başlangıç ihlali)", ('', '1', '', '2x', D(0), NONE, 0, 2), X**2 + 1, 'u(0)')]:
    v = verdict(*args[:6], cand, *args[6:])
    check(f'verify_exact reddeder: {name}', v['status'] == 'rejected' and v['reason'].startswith(start), v['reason'])
v = verdict('1', '', '1', '0', B(0, 1), NONE, sp.sin(X))
check('verify_exact: doğru IVP adayı sembolik doğrulanır', v['status'] == 'verified_symbolic')
o = fem_sheet.exact_only(dict(alpha='x', beta='', gamma='', f='1', x0=0, x1=1, left=D(0), right=D(0)))
check("KR-1: (xu')'=1, u(0)=u(1)=0 için yanlış 'x − 1' artık analitik olarak gönderilmez", o['exact_js'] is None and 'continuous_error' in o, str(o)[:80])
t = fem_sheet.make_sheet(problem('x', '', '', '1', NONE, D(0), 2, 1))
check('kâğıt bölüm 10: kontrol listesi ve kanıt düzeyi yazılı', 'sembolik doğrulandı' in t['text'] and '\\checkmark' in t['text'] and t['exact_status'] == 'verified_symbolic')

# 1g. Birinci mertebe uyumluluk (madde 3): tam yolda kesin, sayısal yolda ölçekli; JS ile aynı kalıntı
fo = fem_sheet.solve_problem(problem('', '1', '', '2x', D(0), NONE, 2, 1, 0, 2))['first_order']
check("u'=2x [0,2] n=2 P1: r₁ = 0 (kesin), uyumlu", fo['residual'] == 0 and fo['compatible'] is True and fo['yb'] == 0)
fo = fem_sheet.solve_problem(problem('', '1', '', '3x^2', D(0), NONE, 2, 1, 0, 2))['first_order']
check("u'=3x² [0,2] n=2 P1: r₁ = −1 (kesin rasyonel), yᵀb = 1, uyumsuz", fo['residual'] == -1 and fo['yb'] == 1 and fo['compatible'] is False, f"{fo['residual']}, {fo['yb']}")
d = problem('', '1', '2', '0', D(1), NONE, 8, 1, 0, 2)
fo = fem_sheet.solve_problem(d)['first_order']
js = json.loads(subprocess.run(['node', str(HERE / 'tests.js'), '--solve-diag', json.dumps(d)], capture_output=True, text=True).stdout)
check("u'+2u=0 n=8: kalıntı SymPy ↔ JS (göreli 1e-12)", abs(float(fo['residual']) - js['residual']) <= 1e-12 * max(1, abs(js['residual'])), f"{float(fo['residual'])} / {js['residual']}")
sheet = fem_sheet.make_sheet(problem('', '1', '', '3x^2', D(0), NONE, 2, 1, 0, 2))['text']
check('kâğıt: 1. mertebede r₁ ve yᵀb satırı, uyumsuzluk açıklaması', 'r_{1} = ' in sheet and 'tam sistem uyumsuz' in sheet)
check("kâğıt: 2. mertebede r satırı yok", 'r_{1} = ' not in fem_sheet.make_sheet(problem('1', '', '', '1', D(0), D(0), 2, 1))['text'])

# 1h. Öğretici türetmeler (kâğıt): şekil fonksiyonu türetmesi hesapla aynı; ödev kıyası adımları; formülasyona göre ayrım
for p_ in (1, 2, 3, 4):
    for h_ in (sp.Integer(1), sp.Rational(1, 2), sp.Rational(1, 3)):
        Ns_, _ = fem_sheet.lagrange_basis(p_, h_)
        try:
            blocks = fem_sheet.shape_function_derivation(p_, h_, Ns_); ok = True
        except RuntimeError:
            ok = False
        check(f'şekil türetmesi P{p_}, h={h_}: Vandermonde yolu hesaptaki N ile aynı', ok and any('\\mathbf{V}\\mathbf{c}' in b_ for b_ in blocks))
hw = fem_sheet.make_sheet(problem('', '1', '', '2x', D(0), NONE, 2, 1, 0, 2))['text']
for label, needle in [('polinom u^e(s) = c₀ + c₁s', 'u^{e}(s) = c_{0} + c_{1} s'), ('koşul c₀ = a₁⁽ᵉ⁾', 'u^{e}(s_{1}) &= c_{0} = a^{(e)}_{1}'),
                      ('çözüm c₁ = â₂ − â₁', 'c_{1} &= - a^{(e)}_{1} + a^{(e)}_{2}'), ('gruplama â₁(1 − s) + â₂(s)', 'a^{(e)}_{1}\\left(1 - s\\right) + a^{(e)}_{2}\\left(s\\right)'),
                      ('global N⁽²⁾ = 2 − x, x − 1', 'N_{1}^{(2)}(x) = 2 - x,\\quad N_{2}^{(2)}(x) = x - 1'), ('u_h⁽¹⁾ = a₁(1 − x) + a₂x', 'u_h^{(1)}(x) = a_{1}\\left(1 - x\\right) + a_{2}\\left(x\\right)'),
                      ('kalıntı R⁽¹⁾ = a₂ − a₁ − 2x', 'R(u_h, x)\\big|_{\\Omega_{1}} = - a_{1} + a_{2} - 2 x'), ('Dirichlet → φ₁ ağırlık uzayında yok', '\\phi_{1} = N_{1}\\ \\text{çıkarıldı: } u(x_0)\\ \\text{verildi (temel koşul)'),
                      ('i = 2 denklemi −a₁/2 + a₃/2 − 2 = 0 (ödev 2. satır)', '- \\frac{a_{1}}{2} + \\frac{a_{3}}{2} - 2 = 0'),
                      ('i = 3 denklemi −a₂/2 + a₃/2 − 5/3 = 0', '- \\frac{a_{2}}{2} + \\frac{a_{3}}{2} - \\frac{5}{3} = 0'),
                      ('yerine koyma (0)(1 − x) + (2/3)(x) = 2x/3', '\\left(0\\right)\\left(1 - x\\right) + \\left(\\frac{2}{3}\\right)\\left(x\\right) = \\frac{2 x}{3}'),
                      ('δ_ij doğrulaması', 'N_i(s_j) = \\delta_{ij}')]:
    check(f'ödev kâğıdı: {label}', needle in hw, needle)
check('1. mertebede kısmi integrasyon yapay olarak eklenmez', 'kısmi integrasyon gerekmez' in hw and 'v\\,(\\alpha u\')\'' not in hw)
bvp = fem_sheet.make_sheet(problem('1', '', '1', 'x', D(0), D(0), 2, 1))['text']
check('2. mertebe: C⁰ gerekçesi + kısmi integrasyon + satır denklemi', 'klasik anlamda tanımsız' in bvp and '(kısmi integrasyon)' in bvp
      and '\\frac{11 a_{2}}{3}' in bvp and 'R_{1}(x) = a_{1} \\left(1 - 2 x\\right) + 2 a_{2} x - x' in bvp and 'dağılımsal kalıntı' in bvp)
ivp = fem_sheet.make_sheet(problem('1', '', '1', '0', B(0, 1), NONE, 2, 1))['text']
check('IVP: akısı bilinmeyen uçtaki test fonksiyonu gerekçesiyle dışarıda', '\\text{bilinmiyor: ağırlık uzayında }' in ivp and 'Petrov–Galerkin' in ivp)

# 1i. Matematiksel çözüm kâğıdı standardı: başlık mertebeye göre; boyutlar; montaj eşlemesi; sınır terimi gerekçesi; hata ölçüleri
check('1. mertebe: "Weak Formulation" başlığı YOK, ağırlıklı kalıntı başlığı var', 'Weak Formulation' not in hw and 'Galerkin Weighted-Residual Formulation' in hw)
check('2. mertebe: "Weak Formulation (Integration by Parts)" var ve ağırlıklı kalıntıdan sonra', 'Weak Formulation (Integration by Parts)' in bvp
      and bvp.index('Galerkin Weighted-Residual') < bvp.index('Weak Formulation (Integration by Parts)'))
check('kâğıt (Heinrich & Pepper): φ_i = N_i ağırlık, u_h = Σ a_j N_j, yerel a_k⁽ᵉ⁾', '\\phi_i(x) = N_i(x)' in hw and 'u_h(x) = \\sum_{j=1}^{3} a_j\\,N_j(x)' in hw
      and 'a^{(e)}_k = a_{(e-1)p+k}' in hw and 'w_i' not in hw and '\\hat{a}' not in hw)
check('kâğıt: a, K, F boyutları ve satır/sütun anlamı', '\\mathbf{K} \\in \\mathbb{R}^{3\\times3}' in hw and 'satır } i \\leftrightarrow \\phi_i' in hw)
check('ödev: montajda K₂₂ = K⁽¹⁾₂₂ + K⁽²⁾₁₁ = 0 ve F₂ = 2/3 + 4/3', 'K_{22} = K^{(1)}_{22} + K^{(2)}_{11}' in hw and 'F_{l,2} = F^{(1)}_{l,2} + F^{(2)}_{l,1}' in hw)
check('2. mertebe: her uçta sınır terimi gerekçesi', bvp.count('x = 0:\\quad') == 1 and bvp.count('x = 1:\\quad') == 1 and '\\phi_i(0) = 0' in bvp)
neu = fem_sheet.make_sheet(problem('-1', '', '', 'sin(pi*x)', D(0), {'type': 'neumann', 'value': '1'}, 2, 2))['text']
check('Neumann ucu: Q = α u\' bilinen akı olarak yazılır', 'Q_1 = \\alpha(1)\\,u\'(1) = -1' in neu and 'Q_1 = -1\\ \\text{(doğal)}' in neu)
check('karşılaştırma: düğüm, noktasal, L² ve H¹ ayrı', all(k in neu for k in ('(düğüm hatası)', '(noktasal hata', '_{L^2}', '_{H^1}')))
em = fem_verify.error_measures(fem_sheet.solve_problem(problem('', '1', '', '2x', D(0), NONE, 2, 1, 0, 2)), fem_sheet.x**2)
check("ödev hata ölçüleri elle hesaba eşit: düğüm 1/3, |e|_H¹ = √(8/9)", abs(em['nodal'] - 1 / 3) < 1e-12 and abs(em['h1'] - (8 / 9) ** 0.5) < 1e-12, f"{em['nodal']} {em['h1']}")
numeric_sheet = fem_sheet.build_sheet(*fem_sheet.run_pipeline(problem('1', '', '', 'abs(x-0.3)', D(0), D(0), 2, 1), symbolic=False))
exact_sheet = fem_sheet.make_sheet(problem('1', '', '', 'x', D(0), D(0), 2, 1))['text']
check('integral kuralı yazılır: sayısal yolda tanh–sinh, tam yolda SymPy', 'tanh–sinh' in numeric_sheet and 'süre sınırı' in numeric_sheet
      and 'SymPy ile tam' in exact_sheet and 'tanh–sinh' not in exact_sheet)

# 1i-b. Kâğıt sunumu: düğüm sayısı ≠ bilinmeyen sayısı; eleman integrali x'ten s'ye dönüşümüyle (J = 1); hata ölçülerinin anlamı
check('kâğıt: N = np + 1 düğüm ile temel koşuldan sonraki bilinmeyen sayısı ayrı', 'temel koşuldan sonra bilinmeyen: } 2' in hw and 'toplam düğüm = bilinmeyen' not in hw)
check('kâğıt: eleman integrali önce x üzerinde, sonra s = x − x_e ile (J = 1)', 'K^{(e)}_{ij} = \\int_{x_e}^{x_e + h}' in hw and 'J = \\frac{dx}{ds} = 1' in hw
      and '\\int_{1}^{2}(\\cdots)\\,dx = \\int_0^{1}(\\cdots)\\big|_{x = s + 1}\\,ds' in hw)
check('kâğıt: her hata ölçüsünün ne ölçtüğü yazılır', 'düğüm hatası: yalnız düğümlerde' in hw and 'türevin hatası' in hw)
same_sheet = fem_sheet.make_sheet(problem('1', '', '1', '1', D(0), D(0), 3, 1))['text']
check('kâğıt: aynı eleman matrisi sessizce atlanmaz (gerekçe + e = 2, 3 satırı)', 'J = 1 \\Rightarrow s \\text{ cinsinden integrand }' in same_sheet
      and 'e = 2,\\ 3:\\qquad \\mathbf{K}^{(e)} = \\mathbf{K}^{(1)}' in same_sheet)

# 1j. V0.1 kalıcı kıyas problemleri (iç regresyon; arayüze girmez). Beklenen değerler ELLE türetildi (motordan değil).
FL = lambda v: {'type': 'flux', 'value': str(v)}
status = lambda sheet: {it['name']: it['status'] for it in sheet['verification']}
R_ = sp.Rational
# Test A: u' + u = 1, u(0) = 0 → u = 1 − e^{−x}. P1, h = 1/4: K^e = ∫N_i(N_j' + N_j) = [[−1/2 + h/3, 1/2 + h/6], [−1/2 + h/6, 1/2 + h/3]], F^e = [h/2, h/2]
dA = problem('', '1', '1', '1', D(0), NONE, 4, 1)
rA = fem_sheet.solve_problem(dA)
check('A: K^(1) = [[−5/12, 13/24], [−11/24, 7/12]], F^(1) = [1/8, 1/8] (elle)',
      rA['elements'][0]['K'] == sp.Matrix([[R_(-5, 12), R_(13, 24)], [R_(-11, 24), R_(7, 12)]]) and list(rA['elements'][0]['F']) == [R_(1, 8), R_(1, 8)])
sA = fem_sheet.make_sheet(dA)
check('A: kalıntı tüm terimleri içerir R = u_h′ + u_h − 1', 'R(u_h, x)\\big|_{\\Omega_{1}} = a_{1} \\left(- 4 x - 3\\right) + a_{2} \\left(4 x + 4\\right) - 1' in sA['text'])
check('A: analitik 1 − e^{−x} sembolik doğrulandı', sA['exact'] == '1 - exp(-x)' and sA['exact_status'] == 'verified_symbolic', sA['exact'])
check('A: atılan satır kalıntısı ve bağımsız yeniden hesabı, tüm doğrulamalar geçti',
      all(v == 'geçti' for v in status(sA).values()) and '1. mertebe atılan satır (bağımsız)' in status(sA), str(status(sA)))
emA = fem_verify.error_measures(fem_sheet.solve_problem(dA), 1 - sp.exp(-fem_sheet.x))
check('A: hata ölçüleri pozitif ve göreli hata tanımlı', emA['l2'] > 0 and emA['h1'] > 0 and 0 < emA['rel_l2'] < 1, str(emA['l2']))
# Test B: u'' = −2, u(0) = u(1) = 0 → u = x(1 − x). P1: K^e = (1/h)[[1, −1], [−1, 1]] = 4·[[1, −1], [−1, 1]], F_l^e = −∫N_i f = [h, h]·1 = [1/4, 1/4]
dB = problem('1', '', '', '-2', D(0), D(0), 4, 1)
rB = fem_sheet.solve_problem(dB)
check('B: K^(1) = [[4, −4], [−4, 4]] (kısmi integrasyon işareti +), F_l^(1) = [1/4, 1/4]',
      rB['elements'][0]['K'] == sp.Matrix([[4, -4], [-4, 4]]) and list(rB['elements'][0]['F']) == [R_(1, 4), R_(1, 4)])
check('B: P1 düğümlerde tam: a = [0, 3/16, 1/4, 3/16, 0]', rB['a'] == [0, R_(3, 16), R_(1, 4), R_(3, 16), 0], str(rB['a']))
emB = fem_verify.error_measures(rB, fem_sheet.x * (1 - fem_sheet.x))
check('B: düğüm hatası 0 ama iç hata > 0 (L² = h²/√120·… > 0, H¹ > 0)', emB['nodal'] < 1e-15 and emB['l2'] > 1e-3 and emB['h1'] > 1e-2, f"{emB['nodal']} {emB['l2']} {emB['h1']}")
dist = fem_verify.distributional_residual(rB)
check('B: düğüm sıçramaları J = u_h′(x⁺) − u_h′(x⁻) = −2h = −1/2 (elle) ve ∫φ_i R_e + J_i = 0',
      all(abs(j + 0.5) < 1e-25 for j in dist['jumps']) and dist['status'] == 'geçti', str(dist['jumps']))
sB = fem_sheet.make_sheet(dB)
check('B: kâğıt eleman içi / dağılımsal / zayıf kalıntıyı ayırır', all(k in sB['text'] for k in ('eleman içi klasik kalıntı', 'dağılımsal kalıntı', 'Weak Formulation (Integration by Parts)'))
      and all(v == 'geçti' for v in status(sB).values()), str(status(sB)))
# Test C: ((1+x)u')' = 4x + 2, u(0) = 0, (1+x)u'(1) = 4 (AKI koşulu) → u = x². P1, h = 1/4: K^(1) = 16·∫_0^{1/4}(1+x)dx·[[1, −1], [−1, 1]] = 9/2·[[1, −1], [−1, 1]]
dC = problem('1+x', '', '', '4x+2', D(0), FL(4), 4, 1)
rC = fem_sheet.solve_problem(dC)
check('C: değişken katsayılı K^(1) = 9/2·[[1, −1], [−1, 1]] (elle)', rC['elements'][0]['K'] == R_(9, 2) * sp.Matrix([[1, -1], [-1, 1]]))
check('C: akı koşulu Q₁ = 4 (u′(1) = 4 olsaydı Q₁ = α(1)·4 = 8 olurdu)', rC['Q']['right'] == 4
      and fem_sheet.solve_problem(problem('1+x', '', '', '4x+2', D(0), Nm(4), 4, 1))['Q']['right'] == 8)
check("C: aynı fiziksel koşul u′(1) = 2 ile aynı a (akı ↔ türev çevrilmeden, ayrı yoldan)",
      rC['a'] == fem_sheet.solve_problem(problem('1+x', '', '', '4x+2', D(0), Nm(2), 4, 1))['a'])
sC = fem_sheet.make_sheet(dC)
check('C: analitik x² denklem + u(0) + α u′(1) = 4 ile doğrulandı', sC['exact'] == 'x**2' and any(c[0] == 'α u′(1) = 4' and c[1] for c in sC['exact_checks']), str(sC['exact_checks']))
check('C: tüm bağımsız doğrulamalar geçti, L² hızı ≈ 2', all(v == 'geçti' for v in status(sC).values()), str(status(sC)))
check('C: x² adayı u′(1) = 4 türev koşuluyla REDDEDİLİR (akı ≠ türev)',
      fem_verify.verify_exact(fem_sheet.solve_problem(problem('1+x', '', '', '4x+2', D(0), Nm(4), 2, 1)), fem_sheet.x ** 2)['status'] == 'rejected')

# 1k. Tanı testleri: beklenen sınıf ve uyarı davranışı
def outcome(data):
    try:
        return 'çözüldü', fem_sheet.make_sheet(data)
    except fem_sheet.Diagnosis as error:
        return (error.category, error.level), str(error)
check('tanı: çözüm yok (saf Neumann, ∫f ≠ akı farkı) → K1 sürekli', outcome(problem('1', '', '', '1', Nm(0), Nm(0), 4, 1))[0] == ('K1', 'sürekli problem'))
check('tanı: sonsuz çözüm (u″ = 0, N/N) → K2 sürekli', outcome(problem('1', '', '', '0', Nm(0), Nm(0), 4, 1))[0] == ('K2', 'sürekli problem'))
check('tanı: saf Neumann uyumlu (∫f = 0) → K2 (tek değil), uyumsuz → K1', outcome(problem('1', '', '', '2x-1', Nm(0), Nm(0), 4, 2))[0][0] == 'K2')
kind, sheetN = outcome(problem('1', '', '9.8', '1', D(0), D(0), 8, 2))
check('tanı: rezonansa yakın (γ = 9.8 ≈ π²) ama tek → çözülür, çözücü kalıntısı geçti', kind == 'çözüldü' and status(sheetN)['çözücü kalıntısı ve koşul sayısı'] == 'geçti')
kind, sheetD = outcome(problem('x', '', '', '1', NONE, D(0), 4, 2))
check('tanı: dejenere uç (α(0) = 0) → reddedilmez, kâğıtta dejenere uyarısı', kind == 'çözüldü' and 'dejenere' in sheetD['text'])
kind, msg = outcome(problem('x-0.5', '', '', '1', D(0), D(0), 4, 1))
check('tanı: işaret değiştiren α, düğüm x* = 1/2 → AYRIK sistem K1, sürekli hüküm yok, ilgisiz saf Neumann örneği yok',
      kind == ('K1', 'ayrık sistem') and 'AYRIK' in msg and "saf Neumann" not in msg, msg[:90])
check('tanı: çözücü başarısızlığı (IVP kaba ağ) → K3', category(problem('1', '', '1', '0', B(0, 1), NONE, 64, 1, 0, 512), symbolic=False) == 'K3')
sH = fem_sheet.make_sheet(problem('1', '', '400', '1', D(0), D(0), 2, 1))
check('tanı: çözünürlüğü yetersiz Helmholtz (k = 20, kh = 10) → dalga uyarısı, hız iddiası yok',
      status(sH)['dalga çözünürlüğü (kh)'] == 'uyarı' and status(sH)['yakınsama çalışması (ifem, analitik çözüme karşı)'] != 'geçti')
sP = fem_sheet.make_sheet(problem('-0.01', '1', '', '1', D(0), D(0), 4, 1))
check('tanı: taşınım baskın (Pe = 12.5) → Péclet uyarısı; analitik yoksa hata/yakınsama "yapılmadı"',
      status(sP)['Péclet sayısı'] == 'uyarı' and status(sP)['yakınsama çalışması'] == 'yapılmadı' and 'hesaplanmadı' in sP['text'])

# 1l. Katman sınırları
v1 = fem_sheet.validate_problem(fem_sheet.define_problem(dC))
check('Katman 2: iz tablosunda her alan iki ayrıştırıcıyla eşdeğer', all(t['ok'] for t in v1['trace']) and {t['item'] for t in v1['trace']} >= {'α(x)', 'f(x)', 'x₁ koşulu', 'operatör'})
check('Katman 2: eşdeğerlik denetimi gerçek farkı yakalar', not fem_sheet.same_expression(fem_sheet.x, fem_sheet.x + sp.Rational(1, 10**6), 0, 1)[0])
check('Katman 2: ondalık girdi iki ayrıştırıcıda eşdeğer (0.05 ↔ 1/20)', fem_sheet.same_expression(fem_sheet.parse('0.05x', 't'), fem_sheet.independent_parse('0.05x'), 0, 1)[0])
try:
    fem_sheet.validate_problem(fem_sheet.define_problem({**dB, 'n': 4.5})); check('Katman 2: n = 4.5 sessizce yuvarlanmaz', False)
except fem_sheet.UserError:
    check('Katman 2: n = 4.5 sessizce yuvarlanmaz', True)
import ast as _ast
verify_src = (HERE / 'fem_verify.py').read_text(encoding='utf-8')
called = {n.func.attr for n in _ast.walk(_ast.parse(verify_src)) if isinstance(n, _ast.Call) and isinstance(n.func, _ast.Attribute)} | \
         {n.func.id for n in _ast.walk(_ast.parse(verify_src)) if isinstance(n, _ast.Call) and isinstance(n.func, _ast.Name)}
check('Katman 4 motoru çağırmaz (fem_engine, lagrange_basis, integral, formulation_data, solve_reduced yok)',
      not called & {'fem_engine', 'lagrange_basis', 'integral', 'formulation_data', 'solve_reduced', 'element_residual', 'shape_function_construction'})
sheet_tree = next(n for n in _ast.parse((HERE / 'fem_sheet.py').read_text(encoding='utf-8')).body if getattr(n, 'name', '') == 'build_sheet')
sheet_calls = {getattr(n.func, 'attr', getattr(n.func, 'id', '')) for n in _ast.walk(sheet_tree) if isinstance(n, _ast.Call)}
check('Katman 5 hesap yapmaz (build_sheet içinde integrate/solve/diff/simplify/Abs/quad/exact_solution/error_measures çağrısı yok)',
      not sheet_calls & {'integrate', 'solve', 'diff', 'simplify', 'Abs', 'quad', 'lambdify', 'exact_solution', 'error_measures', 'fem_engine', 'verify'},
      str(sorted(sheet_calls & {'integrate', 'solve', 'diff', 'simplify', 'Abs', 'quad', 'lambdify'})))

# 1m. Belgelerdeki çalışılmış örnek motorla aynı olmalı (README ve docs/MATEMATIK.md elle yazılır; bir kez aralık ile sayılar
#     ayrışmıştı). Değerler motordan hesaplanır, belgede aynen geçmesi ve u_h'nin ortak düğümde sürekli olması sınanır.
dW = problem('', '1', '', '2x', D(0), NONE, 2, 1, 0, 2)
rW = fem_sheet.solve_problem(dW)
xW = fem_sheet.x
piecesW = [sp.expand(sum(rW['a'][e + k] * rW['N'][k] for k in range(2)).subs(fem_sheet.s, xW - e * rW['h'])) for e in range(2)]
emW = fem_verify.error_measures(rW, xW ** 2)
check('çalışılmış örnek (0, 2): a = [0, 2/3, 4], u_h = 2x/3 | 10x/3 − 8/3, x = 1\'de sürekli',
      rW['a'] == [0, R_(2, 3), 4] and piecesW == [2 * xW / 3, 10 * xW / 3 - R_(8, 3)] and piecesW[0].subs(xW, 1) == piecesW[1].subs(xW, 1), str(piecesW))
check('çalışılmış örnek: düğüm 1/3, L² = √(4/135), H¹ = √(8/9) (elle)', abs(emW['nodal'] - 1 / 3) < 1e-12
      and abs(emW['l2'] - (4 / 135) ** 0.5) < 1e-12 and abs(emW['h1'] - (8 / 9) ** 0.5) < 1e-12, f"{emW['l2']} {emW['h1']}")
for doc in ('README.md', 'docs/MATEMATIK.md'):
    txt = (HERE / doc).read_text(encoding='utf-8')
    check(f'{doc}: örnek motorla aynı (a, u_h, hatalar) ve eski (0, 1) değerleri yok',
          '0,\\ 2/3,\\ 4' in txt and '2x/3$ on $[0, 1]$' in txt and '10x/3 - 8/3$ on $[1, 2]$' in txt and '1/3$ (at $x = 1$)' in txt
          and '0.172133' in txt and '0.942809' in txt and '1/6' not in txt and '0.030429' not in txt and '1/12' not in txt)
math_doc = (HERE / 'docs/MATEMATIK.md').read_text(encoding='utf-8')
r4 = fem_sheet.solve_problem(problem('1', '', '', '1', D(0), NONE, 4, 1)); r4b = fem_sheet.solve_problem(problem('1', '', '', '1', D(0), D(0), 4, 1))
check('MATEMATIK.md: düğüm sayısı ≠ bilinmeyen sayısı (P1, n = 4: 5 düğüm, 4 / 3 bilinmeyen; motorla aynı)',
      '(unknowns)' not in math_doc and '5 nodes; 4 unknowns with one Dirichlet end, 3 with two' in math_doc
      and r4['dofs'] == 5 and len(r4['free']) == 4 and len(r4b['free']) == 3)

# 2. Ayrıştırıcı tarayıcıdakiyle aynı kuralları uygular
x = fem_sheet.x
for text, want in [('1/2x', x / 2), ('2x^2', 2 * x**2), ('-x^2', -x**2), ('e^-x', sp.exp(-x)), ('(x+1)(x-1)', (x + 1) * (x - 1)), ('0.05', sp.Rational(1, 20))]:
    check(f'parse {text}', sp.simplify(fem_sheet.parse(text, 't') - want) == 0)
for bad in ['__import__("os")', 'foo(x)', 'x;1', '1e-8', '2.5E3']:
    try:
        fem_sheet.parse(bad, 't'); check(f'reddeder {bad}', False)
    except fem_sheet.UserError:
        check(f'reddeder {bad}', True)

# 3. Kâğıtlar: SymPy ↔ JavaScript düğüm değerleri, analitik çözüm, LaTeX derlemesi
cases = [
    ("u''+u=x", problem('1', '0', '1', 'x', D(0), D(0), 4, 2), True),
    ("-u''=sin(pi x)", problem('-1', '0', '0', 'sin(pi*x)', D(0), D(0), 3, 1), True),
    ("u''-4u=0, u'(1)=1", problem('1', '0', '-4', '0', D(0), Nm(1), 2, 3), True),
    ("u'+2u=0 (1. mertebe)", problem('0', '1', '2', '0', D(1), NONE, 4, 1), True),
    ("((1+x)u')'=x", problem('1+x', '0', '0', 'x', D(0), Nm(0), 3, 2), True),
    ("-0.05u''+u'=1", problem('-0.05', '1', '0', '1', D(0), D(0), 8, 1), True),
    ("u''-u=1, koşulsuz (doğal)", problem('1', '0', '-1', '1', NONE, NONE, 2, 2), True),
    ("u''+u=cos(x), n=40 p=4", problem('1', '0', '1', 'cos(x)', D(0), D(1), 40, 4, 0, 2), True),
    ("abs yük (sayısal integral)", problem('1', '0', '0', 'abs(x-0.3)', D(0), D(0), 2, 1), False),
    # Başlangıç / son değer problemi: aynı uçta u ve u'
    ("IVP u''+u=0, u(0)=0, u'(0)=1", problem('1', '0', '1', '0', B(0, 1), NONE, 4, 2), True),
    ("IVP sönümlü, n=10 p=3 (ondalıklı dal)", problem('1', '0.4', '4', '0', B(1, 0), NONE, 10, 3, 0, 5), True),
    ("son değer u''-u=0, u(1)=1, u'(1)=0", problem('1', '0', '-1', '0', NONE, B(1, 0), 3, 2), True),
]
for name, data, expect_exact in cases:
    start = time.time()
    sheet = fem_sheet.make_sheet(data)
    seconds = time.time() - start
    gap = max(abs(a - b) for a, b in zip(sheet['nodal'], javascript_solution(data)))
    # Köşeli yükte tarayıcının kuadratürü yaklaşıktır (arayüz bunu uyarıyla bildirir); kâğıt ise integrali tam alır.
    tolerance = 1e-8 if expect_exact else 1e-5
    check(f'{name}: SymPy = JavaScript', gap < tolerance, f'fark {gap:.1e}, {seconds:.1f} sn, {len(sheet["text"])} karakter')
    text = sheet['text']
    check(f'{name}: kâğıt eksiksiz', all(h in text for h in ['## Verilenler', '. Gösterim', '. Ağ', '. Şekil fonksiyonları', '. Yaklaşık çözüm', '. Galerkin Weighted-Residual', '. Eleman matrisleri', '. Montaj', '. Sınır koşulları ve indirgenmiş', '. Çözüm', '. Yaklaşık çözümün yeniden kurulması'])
          and text.index('u_h(x) = \\sum_{j=1}^{') < text.index('R(u_h, x) =') and 'N_i(s_j) = \\delta_{ij}' in text   # R, u_h'den sonra; şekil türetmesi var
          and 'nan' not in text.lower() and 'None' not in text and '\\mathbf{K}' in text)
    if expect_exact:
        check(f'{name}: analitik çözüm bulundu', sheet['exact'] is not None, sheet['exact'] or '')
    if shutil.which('pdflatex'):
        with tempfile.TemporaryDirectory() as folder:
            (Path(folder) / 's.tex').write_text(sheet['tex'], encoding='utf-8')
            run = subprocess.run(['pdflatex', '-no-shell-escape', '-interaction=nonstopmode', '-halt-on-error', 's.tex'], cwd=folder, capture_output=True)
            check(f'{name}: LaTeX derlenir', run.returncode == 0 and (Path(folder) / 's.pdf').exists())

# 4. −u'' = f, lineer eleman: düğüm değerleri analitik çözümle TAM çakışır (bilinen FEM özelliği)
r = fem_sheet.solve_problem(problem('-1', '0', '0', 'sin(pi*x)', D(0), D(0), 4, 1))
exact = lambda t: sp.sin(sp.pi * t) / sp.pi**2
check('lineer eleman düğümde tam', all(sp.simplify(r['a'][i] - exact(sp.Rational(i, 4))) == 0 for i in range(5)))

# 5. Hatalı girdiler anlaşılır mesajla reddedilir
for name, data in [('tekil (koşulsuz u\'\'=1)', problem('1', '0', '0', '1', NONE, NONE, 2, 1)),
                   # 8'den çok serbest bilinmeyen → ondalıklı çözüm dalı; tekillik orada da Türkçe mesaj vermeli
                   ('tekil, sayısal dal (koşulsuz u\'\'=1, n=4 p=2)', problem('1', '0', '0', '1', NONE, NONE, 4, 2)),
                   ('tekil, sayısal dal (iki uçta u\', n=4 p=2)', problem('1', '0', '0', '1', Nm(0), Nm(1), 4, 2)),
                   ('1. mertebe koşulsuz', problem('0', '1', '1', '0', NONE, NONE, 2, 1)),
                   ('α=β=0', problem('0', '0', '1', 'x', D(0), D(0), 2, 1)),
                   ('x₀ ≥ x₁', problem('1', '0', '0', '1', D(0), D(0), 2, 1, 1, 0)),
                   ('IVP + öbür uçta koşul (aşırı belirli)', problem('1', '0', '1', '0', B(0, 1), D(0), 2, 1)),
                   ('IVP, 1. mertebe', problem('0', '1', '1', '0', B(0, 1), NONE, 2, 1)),
                   ("α(x₀)=0 ucunda u' koşulu (sessizce kaybolurdu)", problem('x', '0', '0', '0', Nm(1), D(1), 4, 1)),
                   ('n tam sayı değil', {**problem('1', '0', '0', '1', D(0), D(0), 2, 1), 'n': 'abc'})]:
    try:
        fem_sheet.make_sheet(data); check(f'reddeder: {name}', False)
    except fem_sheet.UserError as error:
        check(f'reddeder: {name}', True, str(error))

# 5b. Korpus düzeltmelerinin regresyonları (verification/ENGINEERING_NOTES.md, D1–D3; sürekli düzey doğrulaması).
#     Hiçbiri uygulamanın kendi çıktısını beklenen değer olarak kullanmaz: beklenenler elle ya da kapalı biçimden gelir.
# D1: serbest bilinmeyen yok (C46) → eskiden HTTP 500. Kâğıt üretilmeli, doğrulama maddeleri bulunmalı, hiçbiri 'başarısız' olmamalı.
zero = fem_sheet.make_sheet(problem('1', '', '', '1', D(0), D(0), 1, 1))
zst = status(zero)
check('D1: 0 serbest bilinmeyen (P1, n=1, iki uçta u) → kâğıt, a = [0, 0], boş ağırlık uzayı yazılı',
      zero['nodal'] == [0.0, 0.0] and 'ağırlık uzayı boş' in zero['text'], str(zero['nodal']))
check('D1: doğrulama maddeleri var, hiçbiri başarısız değil; ifem maddesi geçti',
      bool(zst) and 'başarısız' not in zst.values() and zst.get('bağımsız yeniden montaj ve çözüm (ifem)') == 'geçti', str(zst))
two = fem_sheet.make_sheet(problem('1', '', '', '0', D(2), D(5), 1, 1))
check('D1: birden çok verilen değer (u(0)=2, u(1)=5, 0 serbest) → a = [2, 5]', two['nodal'] == [2.0, 5.0], str(two['nodal']))
p2 = fem_sheet.solve_problem(problem('1', '', '', '1', D(0), D(0), 1, 2))
check("D1 karşı-sınama: P2, n=1 (1 serbest iç düğüm) → a₂ = u(1/2) = −1/8 (u = x(x−1)/2)", p2['a'][1] == sp.Rational(-1, 8), str(p2['a']))
for bad in ('abc', '2*(', '1/0'):   # boş alan bilinçli olarak 0 sayılır (parse), bozuk sayılmaz
    try:
        fem_sheet.make_sheet(problem('1', '', '', '1', D(bad), D(0), 1, 1)); check(f"D1: bozuk verilen değer '{bad}' reddedilir", False)
    except fem_sheet.UserError as error:
        check(f"D1: bozuk verilen değer '{bad}' reddedilir", True, str(error))
boom = lambda *_: 1 / 0
g1, g2 = fem_verify.guarded('deneme', boom), fem_verify.guarded('deneme', boom, pair=True)
check('D1: istisna atan denetim başarısız + istisna metni (asla geçti değil); pair → (madde, None)',
      g1['status'] == 'başarısız' and 'ZeroDivisionError' in g1['detail'] and g2[0]['status'] == 'başarısız' and g2[1] is None, g1['detail'])
def _diagnose():
    raise fem_verify.Diagnosis('K1', 'sürekli problem', 'deneme')
try:
    fem_verify.guarded('deneme', _diagnose); check('D1: Diagnosis yutulmaz (hat durur)', False)
except fem_verify.Diagnosis:
    check('D1: Diagnosis yutulmaz (hat durur)', True)

# D3: karmaşık matriste rank (svd_c); sanal kısım atılmaz. Beklenen rank'lar elle.
I = sp.I
R = fem_verify.numerical_rank
check('D3: gerçek tam rank 2', R(sp.Matrix([[1, 2], [3, 4]])) == 2)
check('D3: karmaşık tam rank 2 (C40 sınır matrisi [[−2i, 0], [iπ/2, 1]])', R(sp.Matrix([[-2 * I, 0], [I * sp.pi / 2, 1]])) == 2)
check('D3: gerçek tam rank eksikliği → 1', R(sp.Matrix([[1, 2], [2, 4]])) == 1)
check('D3: karmaşık tam rank eksikliği [[1, i], [i, −1]] → 1', R(sp.Matrix([[1, I], [I, -1]])) == 1)
check('D3: neredeyse eksik, σ_min ~ 1e-30 (eşik 1e-40 üstü) → 2', R(sp.Matrix([[1, 1], [1, 1 + sp.Rational(1, 10**30)]])) == 2)
check('D3: σ_min ~ 1e-50 (eşik altı) → 1', R(sp.Matrix([[1, 1], [1, 1 + sp.Rational(1, 10**50)]])) == 1)
check('D3: karmaşık ara biçim, gerçek sonuç: tam sanal sütun ölçeklemesi rank değiştirmez',
      R(sp.Matrix([[I, 0], [2 * I, 1]])) == R(sp.Matrix([[1, 0], [2, 1]])) == 2)
corpus = {c['id']: c['input'] for c in json.loads((HERE / 'verification' / 'corpus' / 'cases.json').read_text(encoding='utf-8'))['cases']}
c40 = fem_sheet.make_sheet(corpus['C40'])
u40 = sp.sympify(c40['exact'] or 'nan', locals={'x': x}) if c40['exact'] else None
target = sp.log(1 + x**2) + sp.atan(x) - sp.pi / 4 - sp.log(2)            # elle: (1+x²)u′ = 2x + 1, u(1) = 0
check('D3: C40 tek çözüm, gerçek analitik çözüm ln(1+x²) + atan x − π/4 − ln 2 doğrulandı',
      u40 is not None and not u40.has(I) and sp.simplify(u40 - target) == 0 and c40['exact_status'].startswith('verified'), c40['exact'] or '')

# Sürekli düzey (verification/continuous_problem.py): uygulamanın sınıflandırıcısını çağırmaz. Beklenenler elle.
sys.path.insert(0, str(HERE / 'verification'))
import continuous_problem as CP   # noqa: E402
cp = {k: CP.analyze(corpus[k]) for k in ('C01', 'C46', 'C64', 'C65', 'C66', 'C67', 'C69', 'C70')}
st3 = lambda k: (cp[k]['existence']['status'], cp[k]['uniqueness']['status'], cp[k]['compatibility']['status'])
check('sürekli: C67 rezonans → NONEXISTENCE, uyumluluk ihlali C = ∫ f sin πx = 1/2, iki yöntem tutarlı',
      st3('C67')[0] == 'NONEXISTENCE_CONFIRMED' and st3('C67')[2] == 'COMPATIBILITY_VIOLATED'
      and 'fark C = 1/2' in cp['C67']['compatibility']['evidence'] and 'tutarlı' in cp['C67']['compatibility'].get('cross_check', ''), str(st3('C67')))
check("sürekli: C64 u''=1, N/N → NONEXISTENCE (C = 1)", st3('C64')[0] == 'NONEXISTENCE_CONFIRMED' and 'fark C = 1' in cp['C64']['compatibility']['evidence'])
check("sürekli: C65 u''=2x−1, N/N → EXISTENCE + NONUNIQUENESS, uyumlu", st3('C65') == ('EXISTENCE_CONFIRMED', 'NONUNIQUENESS_CONFIRMED', 'COMPATIBILITY_SATISFIED'), str(st3('C65')))
check("sürekli: C66 u''=1 koşulsuz → NONEXISTENCE", st3('C66')[0] == 'NONEXISTENCE_CONFIRMED')
check('sürekli: C46 düzenli D/D → EXISTENCE + UNIQUENESS', st3('C46')[:2] == ('EXISTENCE_CONFIRMED', 'UNIQUENESS_CONFIRMED'))
check('sürekli: C01 1. mertebe IVP → Yöntem 1 (varlık + teklik)', st3('C01')[:2] == ('EXISTENCE_CONFIRMED', 'UNIQUENESS_CONFIRMED')
      and 'Yöntem 1' in cp['C01']['existence']['method'])
over = CP.analyze(problem('1', '', '1', '0', B(0, 1), D(0), 2, 1))          # 'u ve u′' + öbür uçta u: aşırı belirli
check('sürekli: aşırı belirli koşul takımı Yöntem 1 kapsamı dışında → INCONCLUSIVE (varlık iddiası yok)',
      over['existence']['status'] == 'INCONCLUSIVE' and 'kapsamında değil' in over['existence']['evidence'], over['existence']['evidence'])
blocked = subprocess.run([sys.executable, '-c', "import sys; sys.modules['skfem'] = None; sys.path.insert(0, 'verification'); "
                          "import continuous_problem; print('ok')"], cwd=HERE, capture_output=True, text=True)
check('sürekli düzey doğrulaması scikit-fem olmadan içe aktarılır (tests.py yalnız requirements.txt ile çalışır)',
      blocked.stdout.strip() == 'ok', blocked.stderr.strip()[-120:])
check('sürekli: C69 (iç noktada α = 0) ve C70 (uçta α = 0) → INCONCLUSIVE, kanıtsız hüküm yok',
      st3('C69')[0] == st3('C70')[0] == 'INCONCLUSIVE' and 'CONFIRMED' not in ' '.join(st3('C69') + st3('C70')), f"{st3('C69')} {st3('C70')}")

# 6. Bağımsız doğrulama paketi (verification/): uygulama ↔ bağımsız mpmath FEM ↔ analitik çözüm, 16 benchmark.
#    Beklenen değerler uygulamanın çıktısından değil analitik çözümden gelir; ayrıca temel kayıtla regresyon.
for mode in ('--quick', '--compare'):
    run = subprocess.run([sys.executable, str(HERE / 'verification' / 'run_verification.py'), mode], capture_output=True, text=True)
    check(f'verification {mode}', run.returncode == 0, (run.stdout.strip().splitlines() or ['?'])[-1])

# 7. Performans düzeltmeleri (CI'da HTTP 504: korpus C48/C50; sembolik yolun sessiz zaman aşımı: C31/C32). Sonucu değil
#    süreyi değiştirdikleri için hem eşdeğerlik (eski yolla aynı sayılar) hem süre sınırı sınanır.
import random as _random
import mpmath as _mp
_random.seed(7)
worst, singular_ok = 0.0, True
for trial in range(40):
    n_, kl, ku = _random.randint(1, 25), _random.randint(0, 4), _random.randint(0, 4)
    rows_ = [[_mp.mpf(_random.uniform(-3, 3)) if -kl <= j - i <= ku else _mp.mpf(0) for j in range(n_)] for i in range(n_)]
    if trial % 5 == 0 and n_ > 2:                                    # sıfır köşegen: pivotlama zorunlu
        for i in range(n_):
            rows_[i][i] = _mp.mpf(0)
    b_ = [_mp.mpf(_random.uniform(-1, 1)) for _ in range(n_)]
    try:
        ref, kref = _mp.lu_solve(_mp.matrix(rows_), _mp.matrix(b_)), _mp.cond(_mp.matrix(rows_))
    except ZeroDivisionError:
        ref = None
    try:
        got, kgot = fem_sheet.band_lu_solve(rows_, b_), fem_sheet.band_condition(rows_)
    except ZeroDivisionError:
        got = None
    if (ref is None) != (got is None):
        singular_ok = False
        continue
    if ref is not None:
        scale = max(abs(ref[i]) for i in range(n_))
        worst = max(worst, float(max(abs(ref[i] - got[i]) for i in range(n_)) / scale / (kref * _mp.eps)),
                    float(abs(kref - kgot) / kref / (kref * _mp.eps)))
check('bant LU / κ₁ ↔ mpmath.lu_solve / cond: 40 rastgele bant matris (pivotlamalı), fark yuvarlama düzeyinde (< κε)',
      singular_ok and worst < 1, f'en büyük fark / (κε) = {worst:.1e}')
singular = [[_mp.mpf(1), _mp.mpf(2)], [_mp.mpf(2), _mp.mpf(4)]]
try:
    fem_sheet.band_condition(singular); raised = False
except ZeroDivisionError:
    raised = True
check('bant LU: tekil matris ZeroDivisionError (mpmath ile aynı sözleşme → K1/K2/K3 sınıflandırması)', raised)
t0 = time.time(); r50 = fem_sheet.solve_problem(problem('-1', '', '', '1', D(0), D(0), 64, 4), symbolic=False); dt = time.time() - t0
xs = [float(r50['x0'] + i * r50['h'] / 4) for i in range(r50['dofs'])]          # düğüm aralığı h/p
check("C50 (−u'' = 1, n = 64, p = 4, 255 bilinmeyen): motor < 20 s ve düğümler x(1−x)/2 ile çakışır",
      dt < 20 and max(abs(float(a) - xv * (1 - xv) / 2) for a, xv in zip(r50['a'], xs)) < 1e-12, f'{dt:.1f} s')
t0 = time.time(); r31 = fem_sheet.solve_problem(problem('1/(1+x)', '', '', '1', D(0), D(0), 6, 1)); dt = time.time() - t0
check("C31 α = 1/(1+x) (K'de log): tam yol yerine ondalıklı çözüm, gerekçe yazılı, < 20 s (önce > 10 dk)",
      dt < 20 and not r31['exact'] and r31['solve_note'] == 'K rasyonel değil', f"{dt:.1f} s, {r31['solve_note']}")
r_rat = fem_sheet.solve_problem(problem('1', '', '1', 'x', D(0), D(0), 2, 1))
check('rasyonel K: tam (kesirli) yol değişmedi (a₂ = −3/44, solve_note yok)', r_rat['exact'] and r_rat['solve_note'] is None and r_rat['a'][1] == sp.Rational(-3, 44))
check("kâğıt: ondalıklı çözümde gerekçe a'nın yanında", '(ondalıklı çözüm: K rasyonel değil)' in fem_sheet.make_sheet(problem('1/(1+x)', '', '', '1', D(0), D(0), 2, 1))['text'])
sym_K = fem_sheet.solve_problem(problem('sqrt(1+x)', '', '', '1', D(0), D(0), 2, 2))['elements'][0]
steps = {st['label']: st for st in sym_K['K_steps']}
check("integral önbelleği: β = 0'da K^(1)_12 ve K^(1)_21 aynı değer ve adımlar, kendi etiketleriyle",
      steps['K^{(1)}_{12}']['value'] == steps['K^{(1)}_{21}']['value'] and steps['K^{(1)}_{21}']['label'] == 'K^{(1)}_{21}'
      and steps['K^{(1)}_{12}']['primitive'] == steps['K^{(1)}_{21}']['primitive'])

sys.exit(1 if failures else 0)
