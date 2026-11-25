#!/usr/bin/env python3
"""fem_verify.py — Katman 4: BAĞIMSIZ DOĞRULAMA.

Motorun (fem_sheet.fem_engine) bildirdiği başarıya güvenmez; çıktısını motorun KULLANMADIĞI yollarla sınar:
  - sürekli problem: genel çözüm + koşul sistemi (dsolve, 'ics' yok) → var / yok / tek değil / belirlenemedi
  - analitik çözüm: aday denklemle VE her koşulla (u, u', akı α u') ayrıca sınanır; geçmeyen aday hiçbir yerde kullanılmaz
  - operatör tutarlılığı (yama testi): V_h'deki bir polinom için K a, KULLANICININ operatörü L(u)'nun doğrudan türevle
    hesaplanıp 24 noktalı Gauss ile integrallenmesiyle karşılaştırılır (terim düşmesi, işaret, değişken katsayı, montaj)
  - bağımsız yeniden montaj ve çözüm: verification/ifem.py (Vandermonde tabanı, Gauss–Legendre, 30 basamak LU)
  - çözücü kalıntısı, geriye hata η ve koşul sayısı κ_∞ (30 basamak)
  - dağılımsal kalıntı: eleman içi klasik kalıntı + iç düğümlerdeki akı sıçramaları J_k = [α u_h']
  - hata ölçüleri (düğüm, noktasal, L², H¹; mutlak ve göreli) ve yakınsama çalışması (ifem ile ardışık ağlar)
  - kararlılık göstergeleri (eleman Péclet sayısı, dalga sayısı × h)
Her madde 'passed' / 'failed' / 'warning' / 'not_run' durumu taşır; yapılamayan doğrulama yapılmış gibi sunulmaz.

Bağımlılık kuralı: bu dosya fem_sheet'ten yalnız temel tipleri, sembolleri ve eşikleri alır; motorun hesap fonksiyonlarını
(fem_engine, lagrange_basis, integral, formulation_data, solve_reduced …) ÇAĞIRMAZ. tests.py bunu kaynakta denetler.
Bu, kod bağımsızlığıdır; matematiksel kanıt değildir (ifem aynı formülasyonu bilinçli olarak yeniden gerçekleştirir).

Bölüm haritası (her bölüm yalnız üsttekileri kullanır):
  1. Ortak araçlar ve bağımsız taban (ifem)
  2. Sürekli problem: sınıflandırma ve analitik çözüm doğrulaması
  3. Ayrık sistem: operatör tutarlılığı, bağımsız çözüm, çözücü kalıntısı
  4. Çözüm kalitesi: dağılımsal kalıntı, hata ölçüleri, yakınsama, kararlılık göstergeleri
  5. verify(): hepsini çalıştırır ve özet maddeleri üretir
"""
import sys
import threading
from pathlib import Path

import mpmath
import sympy as sp

from fem_sheet import x, Diagnosis, COND_WARN, COND_SINGULAR, COMPAT_TOL

# ---------- 1. ORTAK ARAÇLAR VE BAĞIMSIZ TABAN ----------
_dps = mpmath.mp.dps
sys.path.insert(0, str(Path(__file__).resolve().parent / 'verification'))
import ifem                     # noqa: E402  bağımsız FEM (uygulamanın hiçbir dosyasını içe aktarmaz)
mpmath.mp.dps = _dps            # ifem içe aktarılırken mp.dps = 30 yapar; uygulamanın kalanı varsayılan hassasiyette kalsın

GAUSS = ifem.NODES              # 24 noktalı Gauss–Legendre, [-1, 1], 30 basamak
REFERENCE_LIMIT = 129           # bağımsız çözüm (yoğun 30 basamak LU) bu kadar bilinmeyene kadar (≈ 2 s); üstünde 'not_run'
SOLVER_LIMIT = 129              # κ_∞ için ters matris bu boyuta kadar
NONSMOOTH = (sp.Abs, sp.Piecewise, sp.sign, sp.Heaviside, sp.Min, sp.Max)
# Sezgisel eşikler (evrensel değil; testlerle sınandı, README'de belgeli):
OPERATOR_TOL = 1e-10            # yama testi, göreli (Gauss 24 noktanın düzgün integrandda verdiği ~1e-25'in çok üstü)
REFERENCE_TOL = 1e-9            # motor ↔ ifem düğüm değerleri, göreli; κ büyükse κ·1e-14'e genişler
BACKWARD_TOL = 1e-12            # geriye hata η
PECLET_WARN = 1.0               # P1 sabit katsayılı modelde salınımsızlık koşulu Pe_h ≤ 1; p > 1 için h/p ile ölçeklenmiş gösterge
PPW_WARN, PPW_SEVERE = 10, 4    # dalga boyu başına düğüm aralığı (2π / (k h/p)); 10 altı faz hatası, 4 altı çözümsüz
RATE_CONSISTENT, RATE_MATCH = 0.2, 0.3   # ardışık hız farkı / beklenen hızdan sapma (asimptotik bölge ölçütü)


# Metin ayrıntılarında düz sayı (kâğıtta \\text{} içine girer; dec()'in LaTeX'li 10^k biçimi orada kullanılamaz).
num = lambda v: f'{float(v):.6g}'
SUBSCRIPT = lambda k: str(k).translate(str.maketrans('0123456789', '₀₁₂₃₄₅₆₇₈₉'))


def lam(expr):
    return sp.lambdify(x, expr, 'mpmath')


def mp(value):
    return mpmath.mpf(str(sp.N(value, 40)))


# Bağımsız taban: ifem'in Vandermonde tabanı t ∈ [0, 1] üzerinde; ikinci türev C matrisinden (motorun ürün formülünden ayrı).
def basis(p):
    _, N, dN, C = ifem.basis(p)
    d2N = lambda t: [sum(C[m, k] * m * (m - 1) * t ** (m - 2) for m in range(2, p + 1)) for k in range(p + 1)]
    return N, dN, d2N


def ifem_problem(r):
    def condition(side):
        kind, value = r['conditions'][side]
        if kind == 'none':
            return ('none',)
        if kind == 'both':
            return ('both', mp(value), mp(r['slopes'][side]))
        return (kind, mp(value))
    return {'alpha': lam(r['alpha']), 'beta': lam(r['beta']), 'gamma': lam(r['gamma']), 'f': lam(r['f']),
            'x0': mp(r['x0']), 'x1': mp(r['x1']), 'order': r['order'], 'left': condition('left'), 'right': condition('right')}


def item(name, status, detail, **values):
    return {'name': name, 'status': status, 'detail': detail, **values}


# ---------- 2. SÜREKLİ PROBLEM: SINIFLANDIRMA VE ANALİTİK ÇÖZÜM ----------
# Sürekli problemin sınıflandırması (yalnız 1. ve 2. mertebe, m = mertebe). Genel çözüm u = u_p + Σ C_i e_i (dsolve, ics
# YOK) problem tanımının getirdiği q koşula yerleştirilir: M C = r (q × m). rank M = rank[M|r] = m → tek; rank[M|r] >
# rank M → çözüm yok; eşit ama < m → tek değil (Fredholm alternatifi). Dejenere uçta açık VARSAYIM: çözümler kapalı
# aralıkta sınırlıdır (FEM'in C⁰ uzayı); o uçta sonsuza giden temel çözümün katsayısı sıfırlanır. İç dejenere noktada
# sınıflandırma yapılmaz (iki yanda ayrı sabitler gerekir) ve bu açıkça söylenir. Sembolik rank kararı 50 basamaklı
# sayısal rankla karşılaştırılır; uyuşmazlıkta sonuç "sayısal" etiketini alır ve kâğıdı engellemez.
def classify_problem(alpha, beta, gamma, f, x0, x1, conditions, slopes, order, ivp, degeneracy, timeout=10):
    box = {}

    def work():
        try:
            box['res'] = _classify(alpha, beta, gamma, f, x0, x1, conditions, slopes, order, ivp, degeneracy)
        except Exception as error:
            box['res'] = {'kind': 'undetermined', 'reason': 'no closed-form general solution found' if isinstance(error, NotImplementedError)
                          else f'symbolic computation did not complete ({type(error).__name__})'}
    thread = threading.Thread(target=work, daemon=True)
    thread.start()
    thread.join(timeout)
    return box.get('res', {'kind': 'undetermined', 'reason': f'timed out ({timeout} s)'})


# Sembolik rank kararının bağımsız sayısal karşılığı (60 basamak SVD). Modül düzeyinde: doğrudan sınanabilsin diye.
# Genel çözüm karmaşık biçimde gelebilir (ör. log(x ± i)); rank sütun ölçeklemesinden bağımsızdır, bu yüzden karmaşık
# SVD kullanılır (tekil değerler gerçek; Aᴴ örtük) ve sanal kısım atılmaz. Gerçek matriste gerçek SVD (aynı sonuç, daha hızlı).
# Eşik göreli: σ > 10⁻⁴⁰·max(σ_max, 1).
def numerical_rank(mat):
    if mat.rows == 0:
        return 0
    with mpmath.workdps(60):
        A = mpmath.matrix(mat.evalf(60).tolist())
        complex_entries = any(mpmath.im(A[i, j]) != 0 for i in range(A.rows) for j in range(A.cols))
        sv = (mpmath.svd_c if complex_entries else mpmath.svd_r)(A, compute_uv=False)
        top = max(abs(v) for v in sv) if len(sv) else 0
        return sum(1 for v in sv if abs(v) > mpmath.mpf('1e-40') * max(top, 1))


def _classify(alpha, beta, gamma, f, x0, x1, conditions, slopes, order, ivp, degeneracy):
    if order not in (1, 2):
        return {'kind': 'unsupported', 'reason': 'only first and second order are classified'}
    if degeneracy['interior']:
        return {'kind': 'undetermined', 'reason': 'the leading coefficient has an interior zero (x* ≈ '
                + ', '.join(f'{float(pt):.4g}' for pt in degeneracy['interior']) + '); classification would need separate constants on each side and is not done in this version'}
    u = sp.Function('u')
    ode = sp.Eq(sp.diff(alpha * u(x).diff(x), x) + beta * u(x).diff(x) + gamma * u(x), f)
    general = sp.dsolve(ode, u(x)).rhs
    constants = sorted(general.free_symbols - {x}, key=str)
    if general.has(sp.Integral) or len(constants) != order:
        return {'kind': 'undetermined', 'reason': 'no closed-form general solution found'}
    particular = general.subs({c: 0 for c in constants})
    basis_fns = [sp.diff(general, c) for c in constants]
    if sp.simplify(general - particular - sum(c * e for c, e in zip(constants, basis_fns))) != 0:
        return {'kind': 'undetermined', 'reason': 'the general solution is not linear in the constants'}
    equations, assumptions = [], []
    for side, point in (('left', x0), ('right', x1)):
        if point in degeneracy['endpoints']:
            direction = '+' if side == 'left' else '-'
            limits = [sp.limit(e, x, point, direction) for e in basis_fns]
            lp = sp.limit(particular, x, point, direction)
            if lp.is_finite is not True or any(L.is_finite is None for L in limits):
                return {'kind': 'undetermined', 'reason': f'boundedness analysis at the degenerate end x = {sp.latex(point)} could not be performed'}
            unbounded = [c for c, L in zip(constants, limits) if L.is_finite is False]
            for c in unbounded:
                equations.append((f'boundedness x={sp.latex(point)}', c))
            if unbounded:                                   # varsayım yalnız gerçekten bir denklem eklediyse beyan edilir
                assumptions.append(f'x = {sp.latex(point)} degenerate end: the solution is assumed bounded (the unbounded fundamental solution is excluded)')
    zeroed = {expr: 0 for label, expr in equations if expr in constants}
    reduced = general.subs(zeroed)
    for side, point in (('left', x0), ('right', x1)):
        kind, value = conditions[side]
        direction = '+' if side == 'left' else '-'
        at = lambda expr: sp.limit(expr, x, point, direction)
        if kind in ('dirichlet', 'both'):
            equations.append((f'u({sp.latex(point)}) = {sp.latex(value)}', at(reduced) - value))
        if kind == 'neumann':
            equations.append((f"u'({sp.latex(point)}) = {sp.latex(value)}", at(reduced.diff(x)) - value))
        if kind == 'flux':
            equations.append((f"\\alpha u'({sp.latex(point)}) = {sp.latex(value)}", at(alpha * reduced.diff(x)) - value))
        if kind == 'both':
            equations.append((f"u'({sp.latex(point)}) = {sp.latex(slopes[side])}", at(reduced.diff(x)) - slopes[side]))
        if kind == 'none' and order == 2 and not ivp and point not in degeneracy['endpoints']:
            equations.append((f"\\alpha u'({sp.latex(point)}) = 0", at(alpha * reduced.diff(x))))   # formülasyonun doğal koşulu
    M, rhs = sp.linear_eq_to_matrix([e for _, e in equations], constants)
    rank_m, rank_a = M.rank(simplify=True), M.row_join(rhs).rank(simplify=True)
    numeric = (numerical_rank(M), numerical_rank(M.row_join(rhs)))   # bağımsız sayısal kontrol (SVD): sembolik rank kararının karşılığı
    method = 'symbolic' if numeric == (rank_m, rank_a) else 'numeric_50'
    if method != 'symbolic':
        rank_m, rank_a = numeric
    m = len(constants)
    evidence = f'rank M = {rank_m}, rank [M|r] = {rank_a}, m = {m} ({len(equations)} conditions)'
    result = {'method': method, 'evidence': evidence, 'assumptions': assumptions,
              'equations': [label for label, _ in equations], 'general': general}
    if rank_a > rank_m:
        return {**result, 'kind': 'no_solution'}
    if rank_m < m:
        null = M.nullspace()[0]
        return {**result, 'kind': 'non_unique', 'null_function': sp.simplify(sum(v * e for v, e in zip(null, basis_fns)))}
    solved = sp.linsolve([e for _, e in equations], constants)
    values = next(iter(solved))
    return {**result, 'kind': 'unique', 'solution': sp.simplify(general.subs(dict(zip(constants, values))))}


def classify(r):
    return classify_problem(r['alpha'], r['beta'], r['gamma'], r['f'], r['x0'], r['x1'], r['conditions'], r['slopes'],
                            r['order'], r['ivp'], r['degeneracy'])


# Sembolik olarak KANITLANMIŞ "çözüm yok" / "tek değil" hattı durdurur (sayısal kanıt yalnız raporlanır).
def continuous_diagnosis(c):
    if c['kind'] not in ('no_solution', 'non_unique') or c.get('method') != 'symbolic':
        return None
    note = ('; ' + '; '.join(c['assumptions'])) if c['assumptions'] else ''
    if c['kind'] == 'no_solution':
        return Diagnosis('K1', 'continuous', f"the conditions are inconsistent with the equation ({c['evidence']}{note}).")
    return Diagnosis('K2', 'continuous', f"u + c·({sp.sstr(c['null_function'])}) is also a solution ({c['evidence']}{note}).")


# Aday analitik çözümün doğrulanması. Sınananlar: (i) R[u] = (αu')' + βu' + γu − f ≡ 0; (ii) her u koşulu; (iii) her u'
# koşulu; (iv) her AKI koşulu α u' = q (kullanıcının yazdığı biçimde; u' koşuluna çevrilmez); (v) yalnız formülasyonun problem
# tanımına kattığı doğal koşul: 2. mertebe sınır değer probleminde "verilmedi" uç ve α(x_b) ≠ 0 → α u' = 0; dejenere
# "verilmedi" uçta akı yerine sınırlılık; (vi) u kapalı aralıkta sonlu. Kanıt düzeyi: önce sembolik (simplify = 0 ya da
# limit tam değer); olmazsa 50 basamakta |değer| ≤ 1e-30·ölçek → "sayısal" (kanıt değil, sayısal destek).
def verify_exact(r, candidate, timeout=8):
    box = {'checks': []}

    # name = (düz metin, LaTeX): arayüz düz metni, kâğıt LaTeX'i kullanır (pdflatex Unicode simge derlemez).
    def record(name, value, symbolic_zero):
        if symbolic_zero:
            box['checks'].append((name, True, 'symbolic'))
            return True
        try:
            with mpmath.workdps(50):
                v = abs(complex(sp.N(value, 50)))
        except Exception:
            box['checks'].append((name, False, 'not_evaluable'))
            return False
        ok = v <= 1e-30 * box['scale']
        box['checks'].append((name, ok, 'numeric_50' if ok else f'= {num(v)} ≠ 0'))
        return ok

    def zero(expr):
        try:
            return sp.simplify(expr) == 0
        except Exception:
            return False

    def work():
        x0, x1 = r['x0'], r['x1']
        samples = [x0 + (x1 - x0) * (1 - sp.cos(sp.pi * (2 * k + 1) / 34)) / 2 for k in range(17)]   # Çebişev iç noktaları
        try:
            values = [abs(complex(sp.N(candidate.subs(x, pt), 30))) for pt in samples]
        except Exception:
            box['checks'].append((('u finite on the interval', 'u\\ \\text{finite}'), False, 'not_evaluable')); return
        box['scale'] = max([1.0] + values)
        residual = sp.diff(r['alpha'] * candidate.diff(x), x) + r['beta'] * candidate.diff(x) + r['gamma'] * candidate - r['f']
        if zero(residual):
            box['checks'].append((('R[u] ≡ 0', 'R[u] \\equiv 0'), True, 'symbolic'))
        else:
            worst = max(samples, key=lambda pt: abs(complex(sp.N(residual.subs(x, pt), 50))))
            if not record(('R[u] ≡ 0', 'R[u] \\equiv 0'), residual.subs(x, worst), False):
                return
        for side, point in (('left', x0), ('right', x1)):
            kind, value = r['conditions'][side]
            direction = '+' if side == 'left' else '-'
            u_b = sp.limit(candidate, x, point, direction)
            if not u_b.is_finite:
                box['checks'].append(((f'u({point}) finite', f'u({sp.latex(point)})\\ \\text{{finite}}'), False, f'limit = {u_b}')); return
            if kind in ('dirichlet', 'both') and not record((f'u({point}) = {value}', f'u({sp.latex(point)}) = {sp.latex(value)}'), u_b - value, zero(u_b - value)):
                return
            if kind in ('neumann', 'both'):
                target = value if kind == 'neumann' else r['slopes'][side]
                du_b = sp.limit(candidate.diff(x), x, point, direction)
                if not record((f"u′({point}) = {target}", f"u'({sp.latex(point)}) = {sp.latex(target)}"), du_b - target, zero(du_b - target)):
                    return
            if kind == 'flux':
                flux = sp.limit(r['alpha'] * candidate.diff(x), x, point, direction)
                if not record((f"α u′({point}) = {value}", f"\\alpha u'({sp.latex(point)}) = {sp.latex(value)}"), flux - value, zero(flux - value)):
                    return
            if kind == 'none' and r['order'] == 2 and not r['ivp']:
                if r['alpha'].subs(x, point) == 0:
                    box['checks'].append(((f'x = {point} degenerate end: u bounded', f'x = {sp.latex(point)}\\ \\text{{degenerate end: }} u\\ \\text{{bounded}}'), True, 'symbolic'))
                else:
                    flux = sp.limit(r['alpha'] * candidate.diff(x), x, point, direction)
                    if not record((f"α u′({point}) = 0", f"\\alpha u'({sp.latex(point)}) = 0"), flux, zero(flux)):
                        return
        box['done'] = True
    thread = threading.Thread(target=work, daemon=True)
    thread.start()
    thread.join(timeout)
    checks = box['checks']
    if 'done' not in box:
        failed = next((c for c in checks if not c[1]), None)
        return {'status': 'rejected' if failed else 'timeout', 'checks': checks,
                'reason': f'{failed[0][0]}: {failed[2]}' if failed else f'verification did not finish within {timeout} s'}
    status = 'verified_symbolic' if all(c[2] == 'symbolic' for c in checks) else 'verified_numeric'
    return {'status': status, 'checks': checks, 'reason': ''}


# Karmaşık biçimde yazılmış (ör. log(x ± i) içeren) ama gerçek olması gereken aday. Sanal kısım atılmaz: önce Im u ≡ 0
# sembolik olarak gösterilmeye çalışılır (x gerçek); olmazsa aralığın 17 Çebişev noktasında 50 basamakta |Im u| ≤ 1e-40·ölçek
# aranır. İkisi de olmazsa aday kullanılmaz. Gerçek kısım sadeleştirilmiş biçimde döndürülür.
def real_form(candidate, x0, x1):
    expanded = sp.expand_complex(candidate)
    re_part, im_part = sp.simplify(sp.re(expanded)), sp.simplify(sp.im(expanded))
    # x gerçekken arg(g(x) + i·c) = atan2(c, g); sabit c ≠ 0 için atan2(c, g) = sign(c)·π/2 − atan(g/c) (her gerçek g için tam).
    def atan2_real(y, g):
        return sp.sign(y) * sp.pi / 2 - sp.atan(g / y) if (y.is_number and y != 0) else sp.atan2(y, g)
    rewrite = lambda e: sp.simplify(e.replace(sp.arg, lambda z: atan2_real(sp.im(z), sp.re(z))).replace(sp.atan2, atan2_real))
    re_part, im_part = rewrite(re_part), rewrite(im_part)   # aynı tam özdeşlik iki kısma da: Im u ≡ 0 sembolik gösterilebilsin
    if im_part == 0:
        return re_part, ''
    points = [x0 + (x1 - x0) * (1 - sp.cos(sp.pi * (2 * k + 1) / 34)) / 2 for k in range(17)]
    try:
        with mpmath.workdps(50):
            scale = max([1.0] + [abs(complex(sp.N(candidate.subs(x, pt), 50))) for pt in points])
            worst = max(abs(float(sp.N(im_part.subs(x, pt), 50))) for pt in points)
    except Exception:
        return None, 'complex-valued candidate: the imaginary part could not be evaluated'
    if worst <= 1e-40 * scale:
        return re_part, ''
    return None, f'complex-valued candidate: the imaginary part could not be shown to vanish (max |Im u| = {num(worst)})'


# Analitik çözüm: sınıflandırmanın genel çözüm + koşul sisteminden bulduğu TEK çözüm. Aday verify_exact'ten geçmezse
# kullanılmaz (u = None); kâğıt ve tarayıcı nedeni yazar ve hata hesaplanmaz.
def exact_solution(r, c=None):
    c = c or r.get('classification') or classify(r)
    if c['kind'] != 'unique':
        return {'u': None, 'status': 'not_found', 'checks': [], 'reason': c.get('reason', 'the continuous problem does not have a unique solution')}
    candidate = c['solution']
    if candidate.has(sp.Integral) or candidate.free_symbols - {x}:
        return {'u': None, 'status': 'not_found', 'checks': [], 'reason': 'not in closed form'}
    if candidate.has(sp.I):
        candidate, why = real_form(candidate, r['x0'], r['x1'])
        if candidate is None:
            return {'u': None, 'status': 'not_found', 'checks': [], 'reason': why}
    verdict = verify_exact(r, candidate)
    return {'u': candidate if verdict['status'].startswith('verified') else None, 'candidate': candidate, **verdict}


# ---------- 3. AYRIK SİSTEM ----------
def _matrix(M):
    return [[mp(v) for v in row] for row in M.tolist()]


# Eleman integrallerini 24 noktalı Gauss ile hem bütün elemanda hem iki yarıda alır; fark kuadratür hatası tahminidir.
# fn(x, t) bir LİSTE döndürür (ör. her yerel test fonksiyonu için bir değer): integrand noktada bir kez hesaplanır.
def _element_gauss(fn, xe, h):
    def rule(a, length):
        total = None
        for t, w in GAUSS:
            xx = a + length * (t + 1) / 2
            values = fn(xx, (xx - xe) / h)
            total = [w * length / 2 * v for v in values] if total is None else [acc + w * length / 2 * v for acc, v in zip(total, values)]
        return total
    whole, left, right = rule(xe, h), rule(xe, h / 2), rule(xe + h / 2, h / 2)
    halves = [l + r for l, r in zip(left, right)]
    return halves, [abs(w - v) for w, v in zip(whole, halves)]


# Yama testi (operatör eşdeğerliği). u_p(x) = Σ_{m ≤ p} t^m, t = (x − x0)/(x1 − x0), global olarak V_h'dedir; düğüm değerleri
# a_j = u_p(x_j) için Galerkin satırları TAM olarak şunu vermelidir (L(u) = (αu')' + βu' + γu, KULLANICININ operatörü):
#   2. mertebe: Σ_j K_ij a_j = φ_i(x1) α(x1) u_p'(x1) − φ_i(x0) α(x0) u_p'(x0) − ∫ φ_i L(u_p) dx        (kısmi integrasyon)
#   1. mertebe: Σ_j K_ij a_j = ∫ φ_i L(u_p) dx
# Yük: F_l,i = −∫ φ_i f (2. mertebe), +∫ φ_i f (1. mertebe). Sağ taraflar motorun hiçbir parçasını kullanmaz (taban: ifem
# Vandermonde; türev: SymPy, kullanıcının katsayılarından; integral: Gauss 24 + yarılama hata tahmini, 30 basamak).
def operator_consistency(r):
    order, n, p, dofs = r['order'], r['n'], r['degree'], r['dofs']
    t = (x - r['x0']) / (r['x1'] - r['x0'])
    u = sum(t ** m for m in range(p + 1))
    L = (sp.diff(r['alpha'] * u.diff(x), x) if order == 2 else 0) + r['beta'] * u.diff(x) + r['gamma'] * u
    with mpmath.workdps(30):
        X0, H = mp(r['x0']), mp(r['h'])
        Nb, _, _ = basis(p)
        Lf, ff, uf = lam(L), lam(r['f']), lam(u)
        rhs, load, est_K, est_F = [mpmath.mpf(0)] * dofs, [mpmath.mpf(0)] * dofs, mpmath.mpf(0), mpmath.mpf(0)
        sgn = -1 if order == 2 else 1
        for e in range(n):
            xe = X0 + e * H
            vals, errs = _element_gauss(lambda xx, tt: [Nv * Lf(xx) for Nv in Nb(tt)] + [Nv * ff(xx) for Nv in Nb(tt)], xe, H)
            for k in range(p + 1):
                rhs[e * p + k] += sgn * vals[k]
                load[e * p + k] += sgn * vals[p + 1 + k]
            est_K, est_F = est_K + sum(errs[:p + 1]), est_F + sum(errs[p + 1:])
        if order == 2:
            du = lam(u.diff(x))
            a_ = lam(r['alpha'])
            rhs[dofs - 1] += a_(mp(r['x1'])) * du(mp(r['x1']))
            rhs[0] -= a_(X0) * du(X0)
        nodes = [X0 + H * mpmath.mpf(i) / p for i in range(dofs)]
        a = [uf(xx) for xx in nodes]
        K, Fl = _matrix(r['K']), [mp(v) for v in r['F_l']]
        Ka = [sum(K[i][j] * a[j] for j in range(dofs)) for i in range(dofs)]
        scale_K = max([mpmath.mpf(1)] + [sum(abs(K[i][j] * a[j]) for j in range(dofs)) for i in range(dofs)])
        scale_F = max([mpmath.mpf(1)] + [abs(v) for v in load])
        gap_K = max(abs(Ka[i] - rhs[i]) for i in range(dofs)) / scale_K
        gap_F = max(abs(Fl[i] - load[i]) for i in range(dofs)) / scale_F
        tol_K = OPERATOR_TOL + 10 * float(est_K / scale_K)
        tol_F = OPERATOR_TOL + 10 * float(est_F / scale_F)
    ok = gap_K <= tol_K and gap_F <= tol_F
    rough = max(float(est_K / scale_K), float(est_F / scale_F)) > 1e-8
    status = 'passed' if ok and not rough else 'warning' if ok else 'failed'
    detail = (f'uₚ = 1 + t + … + tᵖ (p = {p}): max|Ka − ref|/scale = {float(gap_K):.1e}, max|Fₗ − ref|/scale = {float(gap_F):.1e}'
              + (' (non-smooth coefficient/load: large quadrature estimate, tolerance widened)' if rough else ''))
    return item('operator consistency (patch test)', status, detail, probe=sp.sstr(u), gap_K=float(gap_K), gap_F=float(gap_F))


# Motorun çözümünü, AYRI montajlanmış ve 30 basamakta çözülmüş bağımsız FEM'le (ifem) düğüm düğüm karşılaştırır.
# 1. mertebede atılan satırın kalıntısı bağımsız K ve F ile yeniden hesaplanır ve motorunkiyle karşılaştırılır.
def independent_solution(r, condition):
    if r['dofs'] > REFERENCE_LIMIT:
        return item('independent re-assembly and solve (ifem)', 'not_run', f"{r['dofs']} unknowns > {REFERENCE_LIMIT} (time limit)"), None
    with mpmath.workdps(30):
        try:
            ref = ifem.solve(ifem_problem(r), r['n'], r['degree'])
        except ZeroDivisionError:
            return item('independent re-assembly and solve (ifem)', 'failed', 'the independent system is singular (the engine solved it)'), None
        a = [mp(v) for v in r['a']]
        scale = max([mpmath.mpf(1)] + [abs(v) for v in ref['c']])
        gap = float(max(abs(a[i] - ref['c'][i]) for i in range(r['dofs'])) / scale)
    tol = max(REFERENCE_TOL, (condition or 1) * 1e-14)
    rough = any(e.has(*NONSMOOTH) for e in (r['alpha'], r['beta'], r['gamma'], r['f']))
    status = 'passed' if gap <= tol else 'warning' if rough else 'failed'
    detail = (f'Vandermonde basis, Gauss–Legendre 24, 30-digit LU: max|a − a(ifem)|/max|a(ifem)| = {gap:.1e} (tolerance {tol:.0e})'
              + ('; non-smooth coefficient/load: the Gauss rule is approximate, the difference comes from quadrature' if rough and gap > tol else ''))
    return item('independent re-assembly and solve (ifem)', status, detail, gap=gap), ref


def first_order_check(r, ref):
    fo = r['first_order']
    if fo is None:
        return None
    if ref is None:
        return item('first order: dropped row (independent)', 'not_run', 'no independent solution')
    d = fo['row'] - 1
    with mpmath.workdps(30):
        a = [mp(v) for v in r['a']]
        residual = sum(ref['K'][d, j] * a[j] for j in range(r['dofs'])) - ref['F'][d]
        engine = mp(fo['residual'])
        scale = max([mpmath.mpf(1), abs(residual)] + [abs(ref['F'][i]) for i in range(r['dofs'])])
        gap = float(abs(residual - engine) / scale)
    ok = gap <= 1e-9
    return item('first order: dropped row (independent)', 'passed' if ok else 'failed',
                f'r{SUBSCRIPT(fo["row"])}: engine {num(engine)}, with independent K and F {num(residual)} (relative difference {gap:.1e}); '
                + ('the full system is consistent' if fo['compatible'] else 'the full system is INCONSISTENT: a different row choice gives a different a' if fo['compatible'] is False
                   else 'compatibility could not be determined'), residual=float(residual))


# Çözücü kalıntısı ve geriye hata (Rigal–Gaches, ∞-norm): η = ‖F − K a‖ / (‖K‖‖a‖ + ‖F‖); payda sıfırsa η = 0 (kalıntı da
# sıfırsa) ya da ∞. κ_∞ = ‖K‖‖K⁻¹‖ 30 basamakta. Motorun "çözdüm" demesine güvenilmez; kalıntı yeniden hesaplanır.
def solver_check(r):
    if not r['free']:
        return item('solver residual', 'passed', 'no free unknowns'), 1.0
    with mpmath.workdps(30):
        A = mpmath.matrix(_matrix(r['K_ff']))
        b = mpmath.matrix([mp(v) for v in r['F_f']])
        a = mpmath.matrix([mp(r['a'][i]) for i in r['free']])
        res = mpmath.norm(b - A * a, mpmath.inf)
        denom = mpmath.mnorm(A, mpmath.inf) * mpmath.norm(a, mpmath.inf) + mpmath.norm(b, mpmath.inf)
        eta = float(res / denom) if denom != 0 else (0.0 if res == 0 else float('inf'))
        kappa = None
        if A.rows <= SOLVER_LIMIT:
            try:
                # A**-1 ile aynı aritmetik, bant içinde (ifem.band_inverse): yoğun ters O(n³) CI'da süre sınırını aşıyordu
                kappa = float(mpmath.mnorm(A, mpmath.inf) * mpmath.mnorm(ifem.band_inverse(A.tolist()), mpmath.inf))
            except ZeroDivisionError:
                kappa = float('inf')
    status = ('failed' if not eta <= BACKWARD_TOL or (kappa is not None and kappa > COND_SINGULAR)
              else 'warning' if kappa is not None and kappa > COND_WARN else 'passed')
    detail = (f'‖F − K a‖∞ = {num(res)}, η = {eta:.1e} (threshold {BACKWARD_TOL:.0e})'
              + (f', κ∞ (reduced K) = {kappa:.2e}' if kappa is not None else ', κ not computed (size)')
              + (f' > {COND_WARN:.0e}: ill-conditioned, the last digits are unreliable' if kappa is not None and kappa > COND_WARN else ''))
    return item('solver residual and condition number', status, detail, eta=eta, kappa=kappa), kappa


# ---------- 4. ÇÖZÜM KALİTESİ ----------
# Dağılımsal kalıntı (yalnız 2. mertebe). C⁰ u_h için L(u_h) − f = Σ_e R_e χ_{Ω_e} + Σ_k J_k δ(x − x_k), J_k = α(x_k)[u_h']_{x_k}
# (iç düğümlerde türev sıçraması). Galerkin denklemi, iç test fonksiyonu için ∫ φ_i R_e dx + Σ_k φ_i(x_k) J_k = 0 demektir
# (kısmi integrasyonun eleman eleman tersi). Bu özdeşlik motorun a'sı ve BAĞIMSIZ tabanla sınanır; uçlarda α u_h'(x_b) ile
# verilen akı Q arasındaki fark da yazılır (doğal koşul ayrık çözümde yalnız zayıf biçimde sağlanır).
def distributional_residual(r):
    if r['order'] != 2:
        return None
    n, p, dofs = r['n'], r['degree'], r['dofs']
    with mpmath.workdps(30):
        X0, H = mp(r['x0']), mp(r['h'])
        N, dN, d2N = basis(p)
        a = [mp(v) for v in r['a']]
        al, dal, be, ga, ff = lam(r['alpha']), lam(sp.diff(r['alpha'], x)), lam(r['beta']), lam(r['gamma']), lam(r['f'])
        loc = lambda e: a[e * p:(e + 1) * p + 1]
        uh = lambda e, t: sum(c * v for c, v in zip(loc(e), N(t)))
        duh = lambda e, t: sum(c * v for c, v in zip(loc(e), dN(t))) / H
        d2uh = lambda e, t: sum(c * v for c, v in zip(loc(e), d2N(t))) / H ** 2
        R = lambda e, xx, t: dal(xx) * duh(e, t) + al(xx) * d2uh(e, t) + be(xx) * duh(e, t) + ga(xx) * uh(e, t) - ff(xx)
        jumps = []
        for k in range(1, n):
            xk = X0 + k * H
            jumps.append(al(xk) * (duh(k, 0) - duh(k - 1, 1)))
        Rmax = max(abs(R(e, X0 + e * H + H * (t + 1) / 2, (t + 1) / 2)) for e in range(n) for t, _ in GAUSS)
        local, est = [], mpmath.mpf(0)              # ∫ N_k R_e dx her eleman ve yerel k için bir kez (Gauss + yarılama tahmini)
        for e in range(n):
            vals, errs = _element_gauss(lambda xx, tt, e=e: [Nv * R(e, xx, tt) for Nv in N(tt)], X0 + e * H, H)
            local.append(vals)
            est += sum(errs)
        identity, worst = [], mpmath.mpf(0)
        for i in r['rows']:
            if i in (0, dofs - 1):
                continue                                    # uç satırlar sınır terimi taşır; özdeşlik iç satırlar için
            e0, k0 = divmod(i, p)
            value = local[e0][k0] if k0 else local[e0 - 1][p] + local[e0][0] + jumps[e0 - 1]
            identity.append((i + 1, value))
            worst = max(worst, abs(value))
        scale = max([mpmath.mpf(1)] + [abs(j) for j in jumps] + [Rmax * H])
        ends = {}
        for side, e, t, xb in (('left', 0, 0, X0), ('right', n - 1, 1, mp(r['x1']))):
            if r['Q'].get(side) is not None:
                ends[side] = (float(al(xb) * duh(e, t)), float(mp(r['Q'][side])))
    tol = 1e-9 + 10 * float(est / scale)
    rough = float(est / scale) > 1e-8
    ok = float(worst / scale) <= tol
    return item('distributional residual (element interiors + nodal jumps)', 'failed' if not ok else 'warning' if rough else 'passed',
                f'max|Rₑ| = {num(Rmax)}, max|Jₖ| = {num(max([abs(j) for j in jumps] or [0]))}, '
                f'interior rows: max|∫φᵢRₑ + φᵢJ| = {float(worst):.1e}'
                + (' (non-smooth coefficient/load: large quadrature estimate, tolerance widened)' if rough else ''),
                jumps=[float(j) for j in jumps], interior_max=float(Rmax), identity=[(i, float(v)) for i, v in identity],
                ends=ends)


# Hata ölçüleri, FEM yeniden kurulumuyla (eleman eleman polinom, BAĞIMSIZ tabanla): düğüm hatası, noktasal hata (eleman başına
# 200 örnekte en büyük; alt tahmin), L² ve H¹ yarı-normu (mpmath.quad her elemanda ayrı; u_h eleman sınırında kırılır) ve
# göreli karşılıkları (payda sıfırsa göreli hata tanımsız, None).
def error_measures(r, solution):
    p, n = r['degree'], r['n']
    exact_f, exact_d = lam(solution), lam(sp.diff(solution, x))
    N, dN, _ = basis(p)
    with mpmath.workdps(20):
        X0, H = mp(r['x0']), mp(r['h'])
        a = [mp(v) for v in r['a']]
        nodes = [X0 + H * mpmath.mpf(i) / p for i in range(r['dofs'])]
        nodal = max(abs(a[i] - exact_f(nodes[i])) for i in range(r['dofs']))
        u_max = max(abs(exact_f(v)) for v in nodes)
        l2 = h1 = pointwise = norm_l2 = norm_h1 = mpmath.mpf(0)
        midpoints = []
        for e in range(n):
            xe = X0 + e * H
            piece = lambda t, e=e: sum(a[e * p + k] * v for k, v in enumerate(N(t / H)))
            dpiece = lambda t, e=e: sum(a[e * p + k] * v for k, v in enumerate(dN(t / H))) / H
            l2 += mpmath.quad(lambda t: (exact_f(xe + t) - piece(t)) ** 2, [0, H])
            h1 += mpmath.quad(lambda t: (exact_d(xe + t) - dpiece(t)) ** 2, [0, H])
            norm_l2 += mpmath.quad(lambda t: exact_f(xe + t) ** 2, [0, H])
            norm_h1 += mpmath.quad(lambda t: exact_d(xe + t) ** 2, [0, H])
            pointwise = max([pointwise] + [abs(exact_f(xe + H * k / 200) - piece(H * k / 200)) for k in range(201)])
            if e < 6:                                   # ilk iki düğüm arasının ortası: hiçbir p için düğüm değil
                xm = r['x0'] + e * r['h'] + r['h'] / (2 * p)
                midpoints.append((xm, float(exact_f(mp(xm))), float(piece(H / (2 * p)))))
        rel = lambda num, den: float(num / den) if den > 0 else None
        picks = list(range(r['dofs'])) if r['dofs'] <= 13 else sorted(set(round(i * (r['dofs'] - 1) / 10) for i in range(11)))
        xs = [r['x0'] + sp.Rational(i, p) * r['h'] for i in range(r['dofs'])]
        nodal_table = [(i, xs[i], r['a'][i], solution.subs(x, xs[i]), sp.Abs(sp.sympify(r['a'][i]) - solution.subs(x, xs[i])))
                       for i in picks]                  # tam (kesirli) düğüm tablosu: kâğıt yalnız yazar
        return {'nodal_table': nodal_table, 'nodal': float(nodal), 'l2': float(mpmath.sqrt(l2)), 'h1': float(mpmath.sqrt(h1)), 'pointwise': float(pointwise),
                'midpoints': midpoints, 'rel_nodal': rel(nodal, u_max), 'rel_l2': rel(mpmath.sqrt(l2), mpmath.sqrt(norm_l2)),
                'rel_h1': rel(mpmath.sqrt(h1), mpmath.sqrt(norm_h1)), 'rel_pointwise': rel(pointwise, u_max)}


# Kuramsal hız yalnız varsayımları sağlanan sınıfta iddia edilir: 2. mertebe sınır değer problemi, önde gelen katsayı
# aralıkta sıfırlanmıyor, düzgün çözüm → L²: p + 1, H¹: p (Céa + Aubin–Nitsche). Diğer sınıflarda (1. mertebe, başlangıç
# değer problemi, dejenere/işaret değiştiren katsayı) beklenen hız bu sürümde iddia edilmez; yalnız gözlenen hız yazılır.
def expected_rates(r):
    if r['order'] == 2 and not r['ivp'] and r['degeneracy']['kind'] == 'uniform':
        return {'l2': r['degree'] + 1, 'h1': r['degree'], 'basis': 'second-order boundary value problem, α ≠ 0, smooth solution assumed'}
    return None


# Yakınsama çalışması: aynı problem ifem ile ardışık ağlarda (m, 2m, 4m eleman; toplam bilinmeyen ≤ REFERENCE_LIMIT) çözülür,
# hata DOĞRULANMIŞ analitik çözüme karşı Gauss 24 ile ölçülür. Hız r = log(E_h/E_{h/2})/log 2. Ölçüm tabanındaki hata
# (≤ 1e-12·‖u‖) hıza katılmaz. Asimptotik bölge ölçütü (sezgisel): ardışık iki hız RATE_CONSISTENT içinde tutarlı. İki
# değerin azalması tek başına "hız kanıtlandı" sayılmaz.
def convergence_study(r, solution):
    p = r['degree']
    m = r['n']
    while m > 1 and 4 * m * p + 1 > REFERENCE_LIMIT:
        m //= 2
    levels = [k for k in (m, 2 * m, 4 * m) if k * p + 1 <= REFERENCE_LIMIT]
    if len(levels) < 3:
        return item('convergence study', 'not_run', f'three mesh levels do not fit the {REFERENCE_LIMIT}-unknown limit'), None
    uf, duf = lam(solution), lam(sp.diff(solution, x))
    N, dN, _ = basis(p)
    rows = []
    with mpmath.workdps(30):
        prob = ifem_problem(r)
        unorm = None
        for k in levels:
            try:
                sol = ifem.solve(prob, k, p)
            except ZeroDivisionError:
                return item('convergence study', 'not_run', f'the independent system is singular on the n = {k} mesh'), None
            H = sol['h']
            l2 = h1 = un = mpmath.mpf(0)
            for e in range(k):
                xe = sol['x0'] + e * H
                for t, w in GAUSS:
                    tt = (t + 1) / 2
                    xx, wx = xe + H * tt, w * H / 2
                    uh = sum(sol['c'][e * p + j] * v for j, v in enumerate(N(tt)))
                    duh = sum(sol['c'][e * p + j] * v for j, v in enumerate(dN(tt))) / H
                    l2 += wx * (uf(xx) - uh) ** 2
                    h1 += wx * (duf(xx) - duh) ** 2
                    un += wx * uf(xx) ** 2
            unorm = mpmath.sqrt(un)
            rows.append({'n': k, 'h': float(H), 'l2': float(mpmath.sqrt(l2)), 'h1': float(mpmath.sqrt(h1))})
    floor = 1e-12 * max(1.0, float(unorm))
    for norm in ('l2', 'h1'):
        for j in range(1, len(rows)):
            E0, E1 = rows[j - 1][norm], rows[j][norm]
            rows[j][f'rate_{norm}'] = (float(mpmath.log(E0 / E1) / mpmath.log(2)) if E0 > floor and E1 > floor else None)
    expected = expected_rates(r)
    verdicts = []
    for norm in ('l2', 'h1'):
        r1, r2 = rows[1].get(f'rate_{norm}'), rows[2].get(f'rate_{norm}')
        if r1 is None or r2 is None:
            verdicts.append((norm, 'undefined', 'error at the measurement floor (exact approximation or round-off level)'))
            continue
        consistent = abs(r2 - r1) <= RATE_CONSISTENT
        if not consistent:
            verdicts.append((norm, 'not_established', f'successive rates {r1:.2f}, {r2:.2f}: not in the asymptotic range; no rate is claimed'))
        elif expected is None:
            verdicts.append((norm, 'observed', f'observed rate ≈ {r2:.2f} (no theoretical rate is claimed for this class)'))
        else:
            target = expected[norm]
            match = abs(r2 - target) <= RATE_MATCH
            verdicts.append((norm, 'consistent' if match else 'inconsistent', f'observed ≈ {r2:.2f}, expected {target}' + (
                '' if match else ': the expected rate was not observed (the meshes may be pre-asymptotic, e.g. coarse kh/Pe, or the smoothness assumption fails)')))
    # Uyumsuz hız tek başına "yöntem hatalı" kanıtı değildir (ön-asimptotik bölge); hata kanıtı yama testi ve bağımsız
    # çözümden gelir. Bu yüzden durum 'warning'dır, 'failed' değil.
    status = 'passed' if all(v[1] in ('consistent', 'observed', 'undefined') for v in verdicts) else 'warning'
    detail = '; '.join(f'{"L²" if nm == "l2" else "H¹"}: {text}' for nm, _, text in verdicts)
    return item('convergence study (ifem, against the analytical solution)', status, detail, rows=rows, expected=expected,
                verdicts=verdicts), rows


# Kararlılık göstergeleri. Uyarıdır, ret değildir: büyük Péclet ya da kh sonucu GEÇERSİZ kılmaz, yalnız standart Galerkin'in
# bu çözünürlükte güvenilir olmayabileceğini (salınım, faz/kirlilik hatası) söyler.
#   Péclet (2. mertebe, β ≠ 0):  Pe = max |β| h / (2|α|);  P1 sabit katsayılı modelde Pe ≤ 1 salınımsızlık koşulu, p > 1 için
#                                 Pe/p (h/p ile ölçeklenmiş) yaklaşık gösterge olarak kullanılır.
#   Dalga (γ/α > 0 bölgesi):     k = √max(γ/α), dalga boyu başına düğüm aralığı PPW = 2π/(k h/p); < 10 faz hatası, < 4 çözümsüz.
#   1. mertebe (β ≠ 0):           λ = max |γ/β|, λh/p > 1 büyük adım (standart Galerkin bu sınıfta kararlılık garantisi vermez).
def stability_indicators(r):
    order, p, h = r['order'], r['degree'], float(r['h'])
    pts = [r['x0'] + (r['x1'] - r['x0']) * sp.Rational(i, 400) for i in range(401)]
    val = lambda e: [float(lam(e)(mp(pt))) for pt in pts]
    out = []
    if order == 2:
        A, Bv, G = val(r['alpha']), val(r['beta']), val(r['gamma'])
        if any(b != 0 for b in Bv):
            peclets = [abs(b) * h / (2 * abs(a)) for a, b in zip(A, Bv) if a != 0]
            if len(peclets) < len(A):
                out.append(item('Péclet number', 'warning', 'α vanishes in the interval: Pe is undefined there (convection-dominated region)'))
            if peclets:
                pe = max(peclets)
                status = 'warning' if pe / p > PECLET_WARN else 'passed'
                out.append(item('Péclet number', status, f'Pe = max|β|h/(2|α|) = {pe:.3g}, Pe/p = {pe / p:.3g}'
                                + (' > 1: standard Galerkin may oscillate (convection-dominated, coarse mesh)' if status == 'warning' else ' ≤ 1'),
                                value=pe))
        waves = [g / a for a, g in zip(A, G) if a != 0 and g / a > 0]
        if waves:
            k = max(waves) ** 0.5
            ppw = 2 * 3.141592653589793 / (k * h / p)
            status = 'warning' if ppw < PPW_WARN else 'passed'
            out.append(item('wave resolution (kh)', status, f'k = {k:.3g}, kh = {k * h:.3g}, {ppw:.3g} intervals per wavelength'
                            + (' < 4: the mesh cannot resolve the wave; the result is physically unreliable' if ppw < PPW_SEVERE
                               else ' < 10: phase (pollution) error expected' if ppw < PPW_WARN else ''), value=k * h))
    else:
        Bv, G = val(r['beta']), val(r['gamma'])
        rates = [abs(g / b) for b, g in zip(Bv, G) if b != 0]
        if rates and max(rates) > 0:
            lh = max(rates) * h / p
            out.append(item('first-order step indicator', 'warning' if lh > 1 else 'passed',
                            f'λ = max|γ/β| = {max(rates):.3g}, λh/p = {lh:.3g}'
                            + (' > 1: large step; standard Galerkin gives no stability guarantee for this class' if lh > 1 else ''), value=lh))
    return out


# ---------- 5. HEPSİ ----------
# Katman 4'ün girişi: motor çıktısı r (yapılandırılmış veri). Sembolik olarak kanıtlanmış sürekli çözümsüzlük/teksizlik
# hattı durdurur (Diagnosis). Diğer her madde durum bilgisiyle döner; 'not_run' maddeleri açıkça listelenir.
# İsteğe bağlı her denetim yalıtılmış çalışır: biri istisna atarsa o madde 'failed' olarak, istisnanın türü ve metniyle
# raporlanır (asla 'passed' sayılmaz) ve FEM sonucu ile diğer denetimler yine raporlanır. İstisna yutulmaz, görünür kalır.
# pair=True: denetim (madde, ek veri) çifti döndürür; hata halinde ek veri None olur (sonraki denetim onu "yok" sayar).
def guarded(name, fn, *args, pair=False):
    try:
        return fn(*args)
    except Diagnosis:
        raise
    except Exception as error:
        failed = item(name, 'failed', f'error while running the check: {type(error).__name__}: {error}')
        return (failed, None) if pair else failed


def verify(r):
    c = r.get('classification') or classify(r)
    r['classification'] = c
    stop = continuous_diagnosis(c)
    if stop is not None:
        raise stop
    exact = exact_solution(r, c)
    solution = exact['u']
    solver, kappa = guarded('solver residual and condition number', solver_check, r, pair=True)
    reference, ref = guarded('independent re-assembly and solve (ifem)', independent_solution, r, kappa or r.get('condition'), pair=True)
    V = {'classification': c, 'exact': exact, 'operator': guarded('operator consistency (patch test)', operator_consistency, r),
         'reference': reference, 'solver': solver, 'first_order': guarded('first order: dropped row (independent)', first_order_check, r, ref),
         'distributional': guarded('distributional residual (element interiors + nodal jumps)', distributional_residual, r),
         'stability': guarded('stability indicators', stability_indicators, r)}
    if not isinstance(V['stability'], list):
        V['stability'] = [V['stability']]
    if solution is not None:
        try:
            V['errors'] = error_measures(r, solution)
        except Exception as error:                 # hata ölçüsü hesaplanamadı: kâğıt "hesaplanmadı" yazar, madde başarısız
            V['errors'] = None
            V['errors_failure'] = item('error measures', 'failed', f'error while running the check: {type(error).__name__}: {error}')
        V['convergence'], _ = guarded('convergence study (ifem, against the analytical solution)', convergence_study, r, solution, pair=True)
    else:
        why = exact['reason'] or 'no verified analytical solution'
        V['errors'] = None
        V['convergence'] = item('convergence study', 'not_run', f'no verified analytical solution ({why})')
    exact_item = (item('analytical solution', 'passed', 'the equation and all conditions ' + ('symbolically' if exact['status'] == 'verified_symbolic' else 'numerically (50 digits; not a proof)') + ' verified')
                  if solution is not None else
                  item('analytical solution', 'failed' if exact['status'] == 'rejected' else 'not_run',
                       ('candidate rejected: ' if exact['status'] == 'rejected' else 'not found: ') + exact['reason']))
    V['items'] = [exact_item, V['operator'], V['reference'], V['solver']] + [v for v in (V['first_order'], V['distributional']) if v] \
        + [V['convergence']] + V['stability'] + ([V['errors_failure']] if V.get('errors_failure') else [])
    return V
