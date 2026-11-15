#!/usr/bin/env python3
"""fem_sheet.py — Galerkin sonlu eleman ÇÖZÜM KÂĞIDI üreticisi (SymPy ile tam/sembolik hesap).

Problem:   (α(x) u')' + β(x) u' + γ(x) u = f(x),   x ∈ (x0, x1)
Çıktı:     stdin'den JSON alır, stdout'a JSON yazar: { text, tex, nodal, exact, exact_js, K, fixed }
           (mode = 'exact' ise yalnızca { exact_js }: grafikler için analitik çözüm)
           text = "## başlık / $$formül$$" blokları (sayfada MathJax ile çizilir), tex = derlenebilir LaTeX belgesi.

Neden ayrı bir Python dosyası: integrallerin TEK TEK ve TAM alınması bir bilgisayar cebiri sistemi ister.
Tarayıcıdaki JavaScript çekirdeği yalnızca canlı grafikler için sayısal çözer; kâğıttaki her sayı buradan gelir.

Dosya haritası (yukarıdan aşağıya, her bölüm yalnızca üsttekileri kullanır):
  1. Girdi okuma      — metin → SymPy ifadesi (güvenli, beyaz listeli)
  2. Biçimleyiciler   — sayı / matris / denklem → LaTeX
  3. Hesap            — şekil fonksiyonları, eleman integralleri, montaj, sınır koşulları, çözüm
  4. Kâğıt            — hesap sonuçlarını sabit sırayla bloklara döker
  5. LaTeX belgesi ve komut satırı girişi
"""
import json
import os
import re
import sys
import threading

import mpmath
import sympy as sp

x, s = sp.symbols('x s', real=True)

# ---------- 1. GİRDİ OKUMA ----------
# Yalnızca bu adlara izin verilir; metin eval'den önce simge simge denetlendiği için başka hiçbir şey çalıştırılamaz.
NAMES = {
    'x': x, 'pi': sp.pi, 'e': sp.E, 'sin': sp.sin, 'cos': sp.cos, 'tan': sp.tan, 'asin': sp.asin, 'acos': sp.acos,
    'atan': sp.atan, 'sinh': sp.sinh, 'cosh': sp.cosh, 'tanh': sp.tanh, 'exp': sp.exp, 'ln': sp.log, 'log': sp.log,
    'sqrt': sp.sqrt, 'abs': sp.Abs,
}
DEGREE_NAMES = {1: 'lineer', 2: 'kuadratik', 3: 'kübik', 4: 'kuartik'}
MATRIX_CAP = 9        # bundan büyük matrisler ilk 4 / son 2 satır-sütunla gösterilir
DETAIL_LIMIT = 8      # eleman sayısı bunu aşarsa integraller yalnızca ilk 2 ve son elemanda tek tek yazılır
EXACT_SOLVE_LIMIT = 8 # serbest bilinmeyen sayısı bunu aşarsa sistem ondalıklı çözülür (tam çözüm sayfalarca kesir üretir)


class UserError(ValueError):
    """Kullanıcının düzeltebileceği girdi hatası; mesajı aynen arayüzde gösterilir."""


# Teşhis kategorileri (tarayıcıdaki CATEGORY ile aynı adlar; iki çekirdek bağımsız karar verir). Tek bir genel "tekil"
# mesajı yerine NEDEN ve DÜZEY (sürekli problem / ayrık sistem) söylenir; ayrık tekillik, diferansiyel problemin
# çözümsüz olduğu anlamına gelmez ve öyle sunulmaz.
CATEGORY = {'K1': 'Çözüm yok (uyumsuz sistem)', 'K2': 'Tek olmayan çözüm (boş uzay ≠ {0})',
            'K3': 'Neredeyse tekil / kötü koşullu', 'K4': 'Dejenere katsayı', 'K5': 'Sayısal çözücü hatası'}


class Diagnosis(UserError):
    """Kategorili teşhis; UserError olduğu için mevcut hata yolu (Türkçe mesaj, sunucu 400) aynen çalışır."""
    def __init__(self, category, level, detail):
        super().__init__(f'{CATEGORY[category]} — {level}: {detail}')
        self.category, self.level = category, level


# Sezgisel sayısal toleranslar (evrensel sınır DEĞİL; göreli tanımlı, tarayıcıdakilerle aynı değerler, testlerle sınandı).
RANK_TOL, COMPAT_TOL, COND_WARN, COND_SINGULAR = 1e-12, 1e-10, 1e10, 1e14
SVD_ANALYSIS_LIMIT = 65     # bundan büyük sayısal tekil sistemde mpmath SVD çok yavaş: analiz yapılmadığı açıkça söylenir
IVP_COARSE_DETAIL = ("sürekli problem tek çözümlü (lineer başlangıç değer problemi), ama ayrık sistem tekil ya da çok kötü "
                     "koşullu: ağ bu çözüm için fazla kaba (salınımlı çözümde ωh küçük olmalı; p = 1 için ωh < 2√3). n ya da p'yi artır.")


def parse(text, label):
    """'2x^2 + sin(pi*x)' → SymPy. Tarayıcıdaki ayrıştırıcıyla aynı kurallar: ^ üs, örtük çarpım (2x), ondalıklar tam kesir."""
    text = (str(text) if text is not None else '').strip().lower() or '0'
    if re.search(r'\d\.?e[+-]?\d', text):        # "1e-8" = 1·e − 8 okunurdu (e Euler sayısı): sessizce yanlış katsayı
        raise UserError(f'{label}: bilimsel gösterim (ör. 1e-8) desteklenmiyor; "e" Euler sayısıdır. 10^-8 biçiminde yaz.')
    out, previous_is_operand = [], False
    for token in re.findall(r"\d+\.?\d*|\.\d+|[a-z]+|\S", text):
        is_number = token[0].isdigit() or (token[0] == '.' and len(token) > 1)
        is_name = token.isalpha()
        if is_name and token not in NAMES:
            raise UserError(f'{label}: tanınmayan simge "{token}"')
        if not (is_number or is_name or token in list('+-*/^()')):
            raise UserError(f'{label}: beklenmeyen simge "{token}"')
        if previous_is_operand and (is_number or is_name or token == '('):
            out.append('*')                                   # örtük çarpım: 2x, 3sin(x), (x+1)(x-1)
        out.append(f"Rational('{token}')" if is_number else '**' if token == '^' else token)
        previous_is_operand = is_number or token == ')' or token in ('x', 'pi', 'e')
    try:
        return sp.sympify(eval(''.join(out), {'__builtins__': {}}, {**NAMES, 'Rational': sp.Rational}))
    except Exception:
        raise UserError(f'{label}: ifade okunamadı')


def parse_number(text, label):
    value = parse(text, label)
    if value.free_symbols or not value.is_finite:
        raise UserError(f'{label}: sabit bir sayı olmalı')
    return value


# ---------- 2. BİÇİMLEYİCİLER ----------
def dec(value, digits=6):
    """Ondalık gösterim: 6 anlamlı basamak, çok büyük/küçükte a×10^b."""
    value = float(value)
    if abs(value) < 1e-14:
        return '0'
    mantissa, exponent = f'{value:.{digits - 1}e}'.split('e')
    if -4 <= int(exponent) <= 5:
        return f'{float(f"{value:.{digits}g}"):.{digits}g}'
    return f'{float(mantissa):g}\\times 10^{{{int(exponent)}}}'


def tex(value):
    """Tam değer tam yazılır (kesir, π …); ondalıklı değer ondalık yazılır."""
    value = sp.sympify(value)
    if not value.has(sp.Float):
        return sp.latex(value)
    if value.is_number:
        return dec(value)
    return sp.latex(value.xreplace({number: sp.Float(number, 6) for number in value.atoms(sp.Float)}))   # ondalık katsayılı polinom


def tex_value(value):
    """Bir sonucun tam hali; rasyonel değilse yanına ondalık karşılığı."""
    value = sp.sympify(value)
    if value.has(sp.Float) or value.is_Rational:
        return tex(value)
    return f'{sp.latex(value)} \\approx {dec(value.evalf())}'


def visible(count):
    return list(range(count)) if count <= MATRIX_CAP else [0, 1, 2, 3, None, count - 2, count - 1]


def tex_matrix(rows, fmt=tex):
    R, C = visible(len(rows)), visible(len(rows[0]))

    def cell(r, c):
        if r is None and c is None:
            return '\\ddots'
        if r is None:
            return '\\vdots'
        if c is None:
            return '\\cdots'
        return fmt(rows[r][c])
    return '\\begin{bmatrix} ' + ' \\\\ '.join(' & '.join(cell(r, c) for c in C) for r in R) + ' \\end{bmatrix}'


def tex_column(values, fmt=tex):
    return tex_matrix([[v] for v in values], fmt)


def tex_row(values, fmt=tex):
    return tex_matrix([list(values)], fmt) + '^{\\mathsf{T}}'


def tex_terms(pairs):
    """[(katsayı, 'u\\'\\''), …] → '2u'' - u' + (x+1)u' gibi işaretleri düzgün bir toplam. Sıfır katsayılı terim yazılmaz."""
    parts = []
    for coefficient, symbol in pairs:
        coefficient = sp.sympify(coefficient)
        if coefficient == 0:
            continue
        negative = coefficient.could_extract_minus_sign()
        magnitude = -coefficient if negative else coefficient
        if magnitude == 1 and symbol:
            body = symbol
        else:
            body = sp.latex(magnitude)
            if magnitude.is_Add:
                body = f'\\left({body}\\right)'
            body += ('\\,' + symbol) if symbol else ''
        parts.append(('-' if negative else '+', body))
    if not parts:
        return '0'
    head = ('-' if parts[0][0] == '-' else '') + parts[0][1]
    return head + ''.join(f' {sign} {body}' for sign, body in parts[1:])


def index_range(first, last):
    """'i = 2, …, 2' gibi tek elemanlı aralık yazılmasın: tek indis 'i = 2' olarak yazılır."""
    if last - first < 3:                                    # kısa aralık açık yazılır: "i = 2, 3"
        return 'i = ' + ',\\ '.join(str(i) for i in range(first, last + 1))
    return f'i = {first}, \\dots, {last}'


def tex_equation(alpha, beta, gamma, f):
    """Denklemin kendisi. α sabitse α u'' yazılır; değişkense (α u')' biçimi korunur."""
    if alpha.free_symbols:
        lead = f"\\left({tex_terms([(alpha, chr(39).join(['u', '']))])}\\right)'"
        rest = tex_terms([(beta, "u'"), (gamma, 'u')])
        left = lead if rest == '0' else lead + (' ' if rest.startswith('-') else ' + ') + (rest if not rest.startswith('-') else '- ' + rest[1:])
    else:
        left = tex_terms([(alpha, "u''"), (beta, "u'"), (gamma, 'u')])
    return f'{left} = {sp.latex(f)}'


# ---------- 3. HESAP ----------
def lagrange_basis(degree, h):
    """Yerel koordinat s ∈ [0, h] üzerinde eşit aralıklı düğümlü Lagrange şekil fonksiyonları ve türevleri."""
    nodes = [sp.Rational(k, degree) * h for k in range(degree + 1)]
    N = [sp.expand(sp.prod([(s - nodes[m]) / (nodes[k] - nodes[m]) for m in range(degree + 1) if m != k]))
         for k in range(degree + 1)]
    return N, [sp.diff(Nk, s) for Nk in N]


def product_tex(terms):
    """[(işaret, [çarpanlar])] → '(N_i')(α)(N_j') - (N_i)(γ)(N_j)': integrandın açılmamış, yerine konmuş hali.
    Katsayısı sıfır olan terim ve 1'e eşit çarpan yazılmaz."""
    out = ''
    for sign, factors in terms:
        if any(g == 0 for g in factors):
            continue
        body = ''.join(f'\\left({sp.latex(g)}\\right)' for g in factors if g != 1) or '1'
        out += (' - ' if sign < 0 else ' + ' if out else '') + body
    return out or '0'


def integral(label, terms, h, symbolic):
    """Tek bir integrali alır; (kâğıda yazılacak LaTeX satırı, değer) döndürür.
    Satır düzeni: yerine konmuş çarpım → açılmış integrand → ilkel fonksiyon → değer (polinomda hepsi yazılır).
    SymPy kapalı biçim bulamazsa integral sayısal alınır ve '≈' ile işaretlenir."""
    integrand = sum(sign * sp.prod(factors) for sign, factors in terms)
    upper = sp.latex(h)
    raw = f'{label} &= \\int_{{0}}^{{{upper}}}\\left[{product_tex(terms)}\\right]ds'
    wrap = lambda rest: f'\\begin{{aligned}} {raw} \\\\ &{rest} \\end{{aligned}}'
    polynomial = integrand.is_polynomial(s)
    if polynomial:
        integrand = sp.expand(integrand)
    if integrand == 0:
        return wrap('= 0'), sp.Integer(0)
    if symbolic:
        try:
            if polynomial:
                primitive = sp.integrate(integrand, s)
                value = sp.simplify(primitive.subs(s, h) - primitive.subs(s, 0))
                return wrap(f'= \\int_{{0}}^{{{upper}}}\\left({sp.latex(integrand)}\\right)ds = '
                            f'\\left[{sp.latex(primitive)}\\right]_{{0}}^{{{upper}}} = {tex_value(value)}'), value
            value = sp.simplify(sp.integrate(integrand, (s, 0, h)))
            if not value.has(sp.Integral) and value.is_finite and value.is_real:
                return wrap(f'= {tex_value(value)}'), value
        except Exception:
            pass
    value = sp.Float(mpmath.quad(sp.lambdify(s, integrand, 'mpmath'), [0, float(h)]), 15)
    return wrap(f'\\approx {dec(value)}'), value


# Tam (rasyonel/sembolik) yolda tekil sistemin sınıflandırılması: rank ve uyumluluk KESİN (tolerans yok).
def classify_discrete_exact(K_ff, F_f, ivp, has_dirichlet):
    if ivp:
        return Diagnosis('K3', 'ayrık sistem', IVP_COARSE_DETAIL)
    rank, rank_aug, size = K_ff.rank(simplify=True), K_ff.row_join(F_f).rank(simplify=True), K_ff.shape[1]
    if rank_aug > rank:
        return Diagnosis('K1', 'ayrık sistem', f'rank K = {rank} < rank [K|F] = {rank_aug} (tam aritmetik): sağ taraf uyumsuz. '
                                                "Örn. saf Neumann'da ∫f ≠ akı farkı.")
    null = K_ff.nullspace()[0]
    shape = 'sabit vektör (u + c de çözüm)' if len(set(sp.simplify(v / null[0]) for v in null)) == 1 else 'sabit olmayan bir vektör'
    advice = '' if has_dirichlet else " Bir uçta u değeri vermek ya da γ'yı değiştirmek boş uzayı kaldırabilir."   # öneri yalnız u hiç verilmemişse anlamlı
    return Diagnosis('K2', 'ayrık sistem', f'rank K = rank [K|F] = {rank} < {size} (tam aritmetik): boş uzay boyutu {size - rank}, '
                                            f'örnek: {shape}.{advice}')


# Sayısal yolda tekil/neredeyse tekil sistem: mpmath SVD ile rank (σ ≤ RANK_TOL·σ_max) ve normwise geriye hata
# η = ‖b − Ax̂‖/(‖A‖‖x̂‖ + ‖b‖) ile uyumluluk. Büyük sistemde analiz yapılmadığı açıkça söylenir (iddia yok).
def classify_discrete_numeric(A, b, ivp, condition):
    if ivp:
        return Diagnosis('K3', 'ayrık sistem', IVP_COARSE_DETAIL)
    if A.rows > SVD_ANALYSIS_LIMIT:
        return Diagnosis('K3', 'ayrık sistem', f'sayısal olarak tekil (κ₁ ≈ {condition:.1e}); {A.rows} bilinmeyende rank/boş uzay analizi yapılmadı.')
    U, S, V = mpmath.svd_r(A)
    sig = [S[i] for i in range(len(S))]
    smax = max(sig)
    keep = [v > RANK_TOL * smax for v in sig]
    rank = sum(keep)
    xhat = mpmath.matrix(A.cols, 1)
    for j, k in enumerate(keep):
        if k:
            coef = sum(U[i, j] * b[i] for i in range(A.rows)) / sig[j]
            for i in range(A.cols):
                xhat[i] += coef * V[j, i]
    eta = float(mpmath.norm(b - A * xhat) / (smax * mpmath.norm(xhat) + mpmath.norm(b)))
    if rank < A.cols:
        if eta <= COMPAT_TOL:
            return Diagnosis('K2', 'ayrık sistem', f'rank {rank} < {A.cols}; sağ taraf uyumlu (η = {eta:.1e}), boş uzay boyutu {A.cols - rank}.')
        return Diagnosis('K1', 'ayrık sistem', f'rank {rank} < {A.cols}; sağ taraf değer uzayında değil (η = {eta:.1e} > {COMPAT_TOL}).')
    return Diagnosis('K3', 'ayrık sistem', f'sayısal olarak tekil: tam rank ama κ₁ ≈ {condition:.1e}; sonuç güvenilmez.')


# Önde gelen katsayının (2. mertebede α, 1. mertebede β) sıfırları. Neden: α'nın sıfırı ya da işaret değişimi standart
# H¹ kuramını (düzgün eliptiklik |α| ≥ α₀ > 0) bozar; ama TEK BAŞINA problemin çözümsüz olduğunu göstermez. Bu yüzden sonuç
# yalnız sınıflandırma ve uyarıdır, ret değildir. Sıfırlar mümkünse sembolik (solveset), değilse sayısal bulunur ve
# hangi yolun kullanıldığı kaydedilir. Sayısal yol çok dar bir değme sıfırını kaçırabilir (bilinen sınır).
def coefficient_degeneracy(lead, x0, x1):
    span = x1 - x0
    if not lead.free_symbols:
        return {'kind': 'uniform', 'points': [], 'endpoints': [], 'interior': [], 'method': 'sembolik'}
    points, method = None, 'sayısal'
    try:
        if lead.is_rational_function(x):
            found = sp.solveset(sp.numer(sp.together(lead)), x, sp.Interval(x0, x1))
            if isinstance(found, sp.FiniteSet):
                points, method = sorted(found, key=float), 'sembolik'
    except Exception:
        points = None
    if points is None:
        fn = sp.lambdify(x, lead, 'mpmath')
        xs = [x0 + span * sp.Rational(i, 2000) for i in range(2001)]
        vals = [float(fn(mpmath.mpf(float(t)))) for t in xs]
        top = max(abs(v) for v in vals) or 1.0
        points = []
        for i in range(2000):
            if vals[i] == 0 or abs(vals[i]) <= 1e-12 * top:
                points.append(sp.Float(float(xs[i]), 15))
            elif vals[i] * vals[i + 1] < 0:
                points.append(sp.Float(float(mpmath.findroot(fn, (float(xs[i]), float(xs[i + 1])), solver='bisect')), 15))
            elif 0 < i and abs(vals[i]) < abs(vals[i - 1]) and abs(vals[i]) <= abs(vals[i + 1]):
                tmin = mpmath.findroot(lambda t: mpmath.diff(fn, t), float(xs[i]))
                if float(x0) <= tmin <= float(x1) and abs(float(fn(tmin))) <= 1e-12 * top:
                    points.append(sp.Float(float(tmin), 15))
        if abs(vals[-1]) <= 1e-12 * top:
            points.append(sp.Float(float(x1), 15))
        points = sorted(set(points), key=float)
    ends = [pt for pt in points if abs(float(pt - x0)) <= 1e-12 * float(span) or abs(float(pt - x1)) <= 1e-12 * float(span)]
    inner = [pt for pt in points if pt not in ends]
    sign_change = []
    for pt in inner:
        d = float(span) * 1e-7
        left, right = float(lead.subs(x, float(pt) - d)), float(lead.subs(x, float(pt) + d))
        if left * right < 0:
            sign_change.append(pt)
    kind = 'uniform' if not points else 'sign_changing' if sign_change else 'degenerate'
    return {'kind': kind, 'points': points, 'endpoints': [x0 if abs(float(pt - x0)) < abs(float(pt - x1)) else x1 for pt in ends],
            'interior': inner, 'sign_change': sign_change, 'method': method}


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
            box['res'] = {'kind': 'undetermined', 'reason': 'kapalı biçimde genel çözüm bulunamadı' if isinstance(error, NotImplementedError)
                          else f'sembolik hesap tamamlanamadı ({type(error).__name__})'}
    thread = threading.Thread(target=work, daemon=True)
    thread.start()
    thread.join(timeout)
    return box.get('res', {'kind': 'undetermined', 'reason': f'zaman aşımı ({timeout} s)'})


def _classify(alpha, beta, gamma, f, x0, x1, conditions, slopes, order, ivp, degeneracy):
    if order not in (1, 2):
        return {'kind': 'unsupported', 'reason': 'yalnız 1. ve 2. mertebe sınıflandırılır'}
    if degeneracy['interior']:
        return {'kind': 'undetermined', 'reason': 'önde gelen katsayının iç sıfırı var (x* ≈ '
                + ', '.join(dec(pt, 4) for pt in degeneracy['interior']) + '); sınıflandırma iki yanda ayrı sabit gerektirir, bu sürümde yapılmıyor'}
    u = sp.Function('u')
    ode = sp.Eq(sp.diff(alpha * u(x).diff(x), x) + beta * u(x).diff(x) + gamma * u(x), f)
    general = sp.dsolve(ode, u(x)).rhs
    constants = sorted(general.free_symbols - {x}, key=str)
    if general.has(sp.Integral) or len(constants) != order:
        return {'kind': 'undetermined', 'reason': 'kapalı biçimde genel çözüm bulunamadı'}
    particular = general.subs({c: 0 for c in constants})
    basis = [sp.diff(general, c) for c in constants]
    if sp.simplify(general - particular - sum(c * e for c, e in zip(constants, basis))) != 0:
        return {'kind': 'undetermined', 'reason': 'genel çözüm sabitlerde doğrusal değil'}
    equations, assumptions = [], []
    for side, point in (('left', x0), ('right', x1)):
        if point in degeneracy['endpoints']:
            direction = '+' if side == 'left' else '-'
            limits = [sp.limit(e, x, point, direction) for e in basis]
            lp = sp.limit(particular, x, point, direction)
            if lp.is_finite is not True or any(L.is_finite is None for L in limits):
                return {'kind': 'undetermined', 'reason': f'x = {sp.latex(point)} dejenere ucunda sınırlılık analizi yapılamadı'}
            unbounded = [c for c, L in zip(constants, limits) if L.is_finite is False]
            for c in unbounded:
                equations.append((f'sınırlılık x={sp.latex(point)}', c))
            if unbounded:                                   # varsayım yalnız gerçekten bir denklem eklediyse beyan edilir
                assumptions.append(f'x = {sp.latex(point)} dejenere uç: çözüm sınırlı varsayıldı (sonsuza giden temel çözüm dışlandı)')
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
        if kind == 'both':
            equations.append((f"u'({sp.latex(point)}) = {sp.latex(slopes[side])}", at(reduced.diff(x)) - slopes[side]))
        if kind == 'none' and order == 2 and not ivp and point not in degeneracy['endpoints']:
            equations.append((f"\\alpha u'({sp.latex(point)}) = 0", at(alpha * reduced.diff(x))))   # formülasyonun doğal koşulu
    M, r = sp.linear_eq_to_matrix([e for _, e in equations], constants)
    rank_m, rank_a = M.rank(simplify=True), M.row_join(r).rank(simplify=True)
    # Bağımsız sayısal kontrol (50 basamak, SVD): sembolik rank kararlarının sayısal karşılığı
    with mpmath.workdps(60):
        def num_rank(mat):
            if mat.rows == 0:
                return 0
            sv = mpmath.svd_r(mpmath.matrix(mat.evalf(60).tolist()), compute_uv=False)
            top = max(abs(v) for v in sv) if len(sv) else 0
            return sum(1 for v in sv if abs(v) > mpmath.mpf('1e-40') * max(top, 1))
        numeric = (num_rank(M), num_rank(M.row_join(r)))
    method = 'sembolik' if numeric == (rank_m, rank_a) else 'sayısal (50 basamak)'
    if method != 'sembolik':
        rank_m, rank_a = numeric
    m = len(constants)
    evidence = f'rank M = {rank_m}, rank [M|r] = {rank_a}, m = {m} ({len(equations)} koşul)'
    result = {'method': method, 'evidence': evidence, 'assumptions': assumptions,
              'equations': [label for label, _ in equations], 'general': general}
    if rank_a > rank_m:
        return {**result, 'kind': 'no_solution'}
    if rank_m < m:
        null = M.nullspace()[0]
        return {**result, 'kind': 'non_unique', 'null_function': sp.simplify(sum(v * e for v, e in zip(null, basis)))}
    solved = sp.linsolve([e for _, e in equations], constants)
    values = next(iter(solved))
    return {**result, 'kind': 'unique', 'solution': sp.simplify(general.subs(dict(zip(constants, values))))}


# Birinci mertebe uyumluluk (madde 3). Ders kitabı kuralı korunur: u verilen düğümün (d) test fonksiyonu V_h,0'da yoktur,
# denklemi atılır. Tam sistem A a = b (A = K[:, serbest], N × (N−1)) ancak b, A'nın değer uzayındaysa uyumludur; yani her
# yᵀA = 0 için yᵀb = 0. Atılan satır kalıntısı r_d = Σ K_dj a_j − F_d. A₂ (kalan satırlar) tekil değilse sol sıfır uzayı tek
# boyutludur, y = (1 [d'de], −A₂^{-T} a_d) ve yᵀb = −r_d: r_d = 0 ⇔ uyumlu. Tam yolda y = nullspace(Aᵀ) kesin hesaplanır;
# sayısal yolda |yᵀb|/(‖y‖‖b‖) ≤ COMPAT_TOL (sezgisel). Uyumsuzsa başka satır seçimi (başka test uzayı) başka a verir.
def first_order_compatibility(K, F, a, fixed, free, exact):
    (d, u0), = fixed.items()
    A = K.extract(list(range(K.rows)), free)
    b = sp.Matrix([F[i] - K[i, d] * u0 for i in range(K.rows)])
    residual = sum(K[d, j] * a[j] for j in range(K.cols)) - F[d]
    if exact:
        null = (A.T).nullspace()
        y = null[0] / null[0][d] if null and null[0][d] != 0 else (null[0] if null else None)
        yb = sp.simplify((y.T * b)[0]) if y is not None else None
        residual = sp.simplify(residual)
        compatible = residual == 0 if y is not None and len(null) == 1 else None
        scaled = None
    else:
        rows2 = [i for i in range(K.rows) if i != d]
        A2T = mpmath.matrix(A.extract(rows2, list(range(A.cols))).T.evalf(30).tolist())
        z = mpmath.lu_solve(A2T, mpmath.matrix([float(v) for v in A.row(d)]))
        y = [mpmath.mpf(1) if i == d else -z[rows2.index(i)] for i in range(K.rows)]
        bn = [mpmath.mpf(float(v)) for v in b]
        yb = sum(yi * bi for yi, bi in zip(y, bn))
        scaled = float(abs(yb) / (mpmath.sqrt(sum(v * v for v in y)) * mpmath.sqrt(sum(v * v for v in bn))))
        compatible = scaled <= COMPAT_TOL
    return {'row': d + 1, 'residual': residual, 'yb': yb, 'scaled': scaled, 'compatible': compatible, 'exact': exact}


# İndirgenmiş sistemin çözümü (önceki gövde aynen; yalnız ayrı fonksiyona alındı ki ayrık teşhis hatası sınıflandırma
# bilgisiyle birlikte taşınabilsin). Tam yolda kesirli LU, sayısal yolda mpmath LU + κ₁.
def solve_reduced(K_ff, F_f, exact, ivp_unique, has_dirichlet):
    if exact:
        if sp.simplify(K_ff.det()) == 0:
            raise classify_discrete_exact(K_ff, F_f, ivp_unique, has_dirichlet)
        return [sp.simplify(v) for v in K_ff.LUsolve(F_f)], None
    A = mpmath.matrix(K_ff.evalf(30).tolist())
    bvec = mpmath.matrix(F_f.evalf(30).tolist())
    try:                                    # tam tekil matriste mpmath.cond ZeroDivisionError atar: o da sınıflandırılır
        condition = float(mpmath.cond(A)) if K_ff.shape[0] else 1.0
    except ZeroDivisionError:
        condition = float('inf')
    if not condition <= COND_SINGULAR:
        raise classify_discrete_numeric(A, bvec, ivp_unique, condition)
    return [sp.Float(v, 15) for v in mpmath.lu_solve(A, bvec)], condition


def solve_problem(data, symbolic=True, classify=False):
    """Bütün hesabı yapar ve kâğıdın ihtiyaç duyduğu her şeyi tek sözlükte döndürür."""
    alpha, beta = parse(data.get('alpha'), 'α(x)'), parse(data.get('beta'), 'β(x)')
    gamma, f = parse(data.get('gamma'), 'γ(x)'), parse(data.get('f'), 'f(x)')
    x0, x1 = parse_number(data.get('x0'), 'x₀'), parse_number(data.get('x1'), 'x₁')
    try:
        n, degree = int(data.get('n', 4)), int(data.get('degree', 1))
    except (TypeError, ValueError):
        raise UserError('n ve p tam sayı olmalı.')
    if not x0 < x1:
        raise UserError('Aralıkta x₀, x₁\'den küçük olmalı.')
    if not (1 <= n <= 64 and 1 <= degree <= 4):
        raise UserError('n 1–64, p 1–4 aralığında olmalı.')
    order = 2 if alpha != 0 else 1 if beta != 0 else 0
    if order == 0:
        raise UserError('α ve β birlikte sıfır olamaz: bu bir diferansiyel denklem değil.')

    # Koşul türleri: none · dirichlet (u) · neumann (u') · both (aynı uçta u ve u': başlangıç/son değer problemi).
    # 'both' için u' değeri ayrı alanda (slope) gelir; conditions[side] = (tür, u değeri) biçimi korunur ki kâğıt kodu değişmesin.
    conditions, slopes = {}, {}
    for side in ('left', 'right'):
        raw = data.get(side) or {}
        kind = raw.get('type', 'none')
        if kind not in ('none', 'dirichlet', 'neumann', 'both'):
            raise UserError('Bilinmeyen koşul türü.')
        conditions[side] = (kind, parse_number(raw.get('value'), 'koşul değeri') if kind != 'none' else None)
        if kind == 'both':
            slopes[side] = parse_number(raw.get('slope'), "u' koşul değeri")
    kinds = [conditions['left'][0], conditions['right'][0]]
    if order == 1 and ('neumann' in kinds or 'both' in kinds or kinds.count('dirichlet') != 1):
        raise UserError('Birinci mertebe denklem tam olarak bir koşul alır: uçlardan birinde u değeri ver, diğerini boş bırak.')
    # Başlangıç değer problemi: iki koşul aynı uçta. Öbür uca koşul verilirse problem aşırı belirlenir (3 koşul, 2. mertebe).
    ivp = 'both' in kinds
    if ivp and sorted(kinds) != ['both', 'none']:
        raise UserError("Bir uçta hem u hem u' verildiyse (başlangıç değer problemi) diğer uca koşul verilmez.")
    # u' koşulu zayıf forma yalnızca akı Q = α u' üzerinden girer. α(uç) = 0 ise Q = 0 olur ve verilen u' değeri
    # sessizce kaybolur (çözüm koşulu sağlamaz). Bu yüzden u' verilen uçta α ≠ 0 şartı aranır.
    for side, point in (('left', x0), ('right', x1)):
        if conditions[side][0] in ('neumann', 'both') and alpha.subs(x, point) == 0:
            raise UserError("α bu uçta sıfır: u' koşulu denkleme girmez (α u' akısı sıfır olur).")

    # Önde gelen katsayının sıfırları (uyarı düzeyi) ve istenirse sürekli sınıflandırma. Yalnız SEMBOLİK olarak kanıtlanan
    # "çözüm yok" / "tek değil" kâğıdı engeller; sayısal kanıt ya da belirlenemeyen durum yalnız raporlanır.
    degeneracy = coefficient_degeneracy(alpha if order == 2 else beta, x0, x1)
    classification = None
    if classify:
        classification = classify_problem(alpha, beta, gamma, f, x0, x1, conditions, slopes, order, ivp, degeneracy)
        if classification['kind'] in ('no_solution', 'non_unique') and classification['method'] == 'sembolik':
            note = ('; ' + '; '.join(classification['assumptions'])) if classification['assumptions'] else ''
            if classification['kind'] == 'no_solution':
                raise Diagnosis('K1', 'sürekli problem', f"koşullar denklemle uyumsuz ({classification['evidence']}{note}).")
            raise Diagnosis('K2', 'sürekli problem', f"u + c·({sp.sstr(classification['null_function'])}) de çözüm ({classification['evidence']}{note}).")
    ivp_unique = ivp and degeneracy['kind'] == 'uniform'     # varlık-teklik teoremi α ≠ 0 iken geçerli

    h = (x1 - x0) / n
    dofs = n * degree + 1
    N, dN = lagrange_basis(degree, h)
    size = degree + 1
    uniform_K = not (alpha.free_symbols or beta.free_symbols or gamma.free_symbols)
    uniform_F = not f.free_symbols
    # Birinci mertebede (α = 0) kısmi integrasyon yoktur: ∫ N_i (β u' + γ u) dx = ∫ N_i f dx doğrudan yazılır.
    # İkinci mertebe formülün (−1) katıdır; çözüm aynıdır ama kâğıt ders kitabındaki işaretle ve doğru formülasyonla yazılır.
    sign = 1 if order == 2 else -1
    detailed = list(range(n)) if n <= DETAIL_LIMIT else [0, 1, n - 1]

    elements = []                                   # eleman başına: K^e, F_l^e ve (varsa) tek tek integral satırları
    for e in range(n):
        xe = x0 + e * h
        local = lambda g: g.subs(x, xe + s)         # katsayıyı elemanın yerel koordinatına taşı: x = x_e + s
        show = symbolic and e in detailed
        record = {'e': e, 'xe': xe, 'K_lines': [], 'F_lines': []}
        if uniform_K and e > 0:
            record['K'] = elements[0]['K']          # sabit katsayı + eşit eleman boyu → matris her elemanda aynı
        else:
            record['K'] = sp.zeros(size, size)
            for i in range(size):
                for j in range(size):
                    terms = ([(1, [dN[i], local(alpha), dN[j]]), (-1, [N[i], local(beta), dN[j]]), (-1, [N[i], local(gamma), N[j]])]
                             if order == 2 else [(1, [N[i], local(beta), dN[j]]), (1, [N[i], local(gamma), N[j]])])
                    line, value = integral(f'K^{{({e + 1})}}_{{{i + 1}{j + 1}}}', terms, h, show)
                    record['K'][i, j] = value
                    record['K_lines'].append(line)
        if uniform_F and e > 0:
            record['F'] = elements[0]['F']
        else:
            record['F'] = sp.zeros(size, 1)
            for i in range(size):
                line, value = integral(f'F^{{({e + 1})}}_{{l,{i + 1}}}', [(-sign, [N[i], local(f)])], h, show)
                record['F'][i] = value
                record['F_lines'].append(line)
        elements.append(record)

    K, F_l = sp.zeros(dofs, dofs), sp.zeros(dofs, 1)
    for record in elements:                         # montaj: yerel düğüm k → global e·p + k
        base = record['e'] * degree
        for i in range(size):
            F_l[base + i] += record['F'][i]
            for j in range(size):
                K[base + i, base + j] += record['K'][i, j]

    # Sınır terimi F_b = [N_i α u']: yalnızca uç düğümlerde sıfırdan farklıdır. Q0 = (αu')(x0), Q1 = (αu')(x1).
    # Q None = bu uçtaki akı bilinmiyor (u verilmiş uç ya da başlangıç değer probleminde koşulsuz kalan uç).
    # Başlangıç değer probleminde koşulsuz uç "doğal koşul" DEĞİLDİR: orada αu' serbesttir, sıfırlanmaz.
    Q = {}
    if order == 2:
        for side, point in (('left', x0), ('right', x1)):
            kind, value = conditions[side]
            Q[side] = (alpha.subs(x, point) * value if kind == 'neumann' else alpha.subs(x, point) * slopes[side] if kind == 'both'
                       else sp.Integer(0) if kind == 'none' and not ivp else None)
    F_b = sp.zeros(dofs, 1)
    if order == 2:
        F_b[0] -= Q['left'] if Q['left'] is not None else 0
        F_b[dofs - 1] += Q['right'] if Q['right'] is not None else 0

    fixed = {}
    if conditions['left'][0] in ('dirichlet', 'both'):
        fixed[0] = conditions['left'][1]
    if conditions['right'][0] in ('dirichlet', 'both'):
        fixed[dofs - 1] = conditions['right'][1]
    free = [i for i in range(dofs) if i not in fixed]
    # Kullanılan denklemler (test fonksiyonları). Sınır değer probleminde serbest düğümlerin denklemleri (Bubnov–Galerkin).
    # Başlangıç değer probleminde akısı BİLİNMEYEN uçtaki test fonksiyonu çıkarılır, akısı bilinen uçtaki tutulur:
    # bilinmeyen sınır terimi zayıf formdan düşer ve sistem yine kare kalır (n·p denklem, n·p bilinmeyen).
    rows = ([i for i in range(dofs - 1)] if conditions['left'][0] == 'both' else
            [i for i in range(1, dofs)] if conditions['right'][0] == 'both' else free)
    K_ff = K.extract(rows, free)
    F_f = sp.Matrix([F_l[i] + F_b[i] - sum(K[i, j] * v for j, v in fixed.items()) for i in rows])

    # Başlangıç değer probleminde tekillik "çözüm tek değil" demek değildir (sürekli problem tek çözümlüdür): yöntem koşullu
    # kararlıdır; kaba ağda ayrık rezonans (ör. p = 2'de ω²h² = 10) ya da üstel büyüme (p = 1'de ωh > 2√3) sistemi bozar.
    exact = not (K_ff.has(sp.Float) or F_f.has(sp.Float)) and len(free) <= EXACT_SOLVE_LIMIT
    condition = None
    try:
        a_free, condition = solve_reduced(K_ff, F_f, exact, ivp_unique, bool(fixed))
    except Diagnosis as error:              # ayrık tekillik, sürekli sınıflandırmayı silmesin: hatayla birlikte taşınır
        error.classification, error.degeneracy = classification, degeneracy
        raise
    a = [None] * dofs
    for i, v in fixed.items():
        a[i] = v
    for i, v in zip(free, a_free):
        a[i] = v

    reactions = {}                                  # akısı bilinmeyen uçta akı: sistemde kullanılmayan o düğümün denkleminden
    if order == 2:
        if Q['left'] is None:
            reactions['left'] = sp.simplify(F_l[0] - sum(K[0, j] * a[j] for j in range(dofs)))
        if Q['right'] is None:
            reactions['right'] = sp.simplify(sum(K[dofs - 1, j] * a[j] for j in range(dofs)) - F_l[dofs - 1])

    first_order = first_order_compatibility(K, F_l, a, fixed, free, exact) if order == 1 else None

    return dict(alpha=alpha, beta=beta, gamma=gamma, f=f, x0=x0, x1=x1, n=n, degree=degree, order=order, h=h, dofs=dofs, first_order=first_order,
                conditions=conditions, slopes=slopes, ivp=ivp, N=N, dN=dN, elements=elements, detailed=detailed, uniform_K=uniform_K,
                uniform_F=uniform_F, K=K, F_l=F_l, F_b=F_b, Q=Q, fixed=fixed, free=free, rows=rows, K_ff=K_ff, F_f=F_f, a=a,
                exact=exact, reactions=reactions, symbolic=symbolic, condition=condition,
                degeneracy=degeneracy, classification=classification)


# Aday analitik çözümün doğrulanması (madde 1). Sınananlar: (i) R[u] = (αu')' + βu' + γu − f ≡ 0; (ii) her u koşulu;
# (iii) her u' koşulu; (iv) yalnız formülasyonun problem tanımına kattığı akı koşulu: 2. mertebe sınır değer probleminde
# "verilmedi" uç ve α(x_b) ≠ 0 → α u' = 0 (IVP'nin serbest ucunda ve 1. mertebede akı SINANMAZ); dejenere "verilmedi"
# uçta akı yerine sınırlılık (sonlu limit); (v) u kapalı aralıkta sonlu. Kanıt düzeyi: önce sembolik (simplify = 0 ya da
# limit tam değer); olmazsa 50 basamakta |değer| ≤ 1e-30·ölçek → "sayısal" (kanıt değil, sayısal destek).
def verify_exact(r, candidate, timeout=8):
    box = {'checks': []}

    # name = (düz metin, LaTeX): arayüz düz metni, kâğıt LaTeX'i kullanır (pdflatex Unicode simge derlemez).
    def record(name, value, symbolic_zero):
        if symbolic_zero:
            box['checks'].append((name, True, 'sembolik'))
            return True
        try:
            with mpmath.workdps(50):
                v = abs(complex(sp.N(value, 50)))
        except Exception:
            box['checks'].append((name, False, 'değerlendirilemedi'))
            return False
        ok = v <= 1e-30 * box['scale']
        box['checks'].append((name, ok, 'sayısal (50 basamak)' if ok else f'= {dec(v)} ≠ 0'))
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
            box['checks'].append((('u aralıkta sonlu', 'u\\ \\text{sonlu}'), False, 'değerlendirilemedi')); return
        box['scale'] = max([1.0] + values)
        residual = sp.diff(r['alpha'] * candidate.diff(x), x) + r['beta'] * candidate.diff(x) + r['gamma'] * candidate - r['f']
        if zero(residual):
            box['checks'].append((('R[u] ≡ 0', 'R[u] \\equiv 0'), True, 'sembolik'))
        else:
            worst = max(samples, key=lambda pt: abs(complex(sp.N(residual.subs(x, pt), 50))))
            if not record(('R[u] ≡ 0', 'R[u] \\equiv 0'), residual.subs(x, worst), False):
                return
        for side, point in (('left', x0), ('right', x1)):
            kind, value = r['conditions'][side]
            direction = '+' if side == 'left' else '-'
            u_b = sp.limit(candidate, x, point, direction)
            if not u_b.is_finite:
                box['checks'].append(((f'u({point}) sonlu', f'u({sp.latex(point)})\\ \\text{{sonlu}}'), False, f'limit = {u_b}')); return
            if kind in ('dirichlet', 'both') and not record((f'u({point}) = {value}', f'u({sp.latex(point)}) = {sp.latex(value)}'), u_b - value, zero(u_b - value)):
                return
            if kind in ('neumann', 'both'):
                target = value if kind == 'neumann' else r['slopes'][side]
                du_b = sp.limit(candidate.diff(x), x, point, direction)
                if not record((f"u′({point}) = {target}", f"u'({sp.latex(point)}) = {sp.latex(target)}"), du_b - target, zero(du_b - target)):
                    return
            if kind == 'none' and r['order'] == 2 and not r['ivp']:
                if r['alpha'].subs(x, point) == 0:
                    box['checks'].append(((f'x = {point} dejenere uç: u sınırlı', f'x = {sp.latex(point)}\\ \\text{{dejenere uç: }} u\\ \\text{{sınırlı}}'), True, 'sembolik'))
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
                'reason': f'{failed[0][0]}: {failed[2]}' if failed else f'doğrulama {timeout} s içinde bitmedi'}
    status = 'verified_symbolic' if all(c[2] == 'sembolik' for c in checks) else 'verified_numeric'
    return {'status': status, 'checks': checks, 'reason': ''}


# Analitik çözüm: sınıflandırmanın genel çözüm + koşul sisteminden bulduğu TEK çözüm (dsolve 'ics' kullanılmaz; o yol
# tekil uçta koşulu sessizce düşürebiliyordu). Aday verify_exact'ten geçmezse kullanılmaz; kâğıt ve tarayıcı nedeni yazar.
def exact_solution(r):
    c = r.get('classification') or classify_problem(r['alpha'], r['beta'], r['gamma'], r['f'], r['x0'], r['x1'], r['conditions'],
                                                   r['slopes'], r['order'], r['ivp'], r['degeneracy'])
    if c['kind'] != 'unique':
        return {'u': None, 'status': 'not_found', 'checks': [], 'reason': c.get('reason', 'sürekli problem tek çözümlü değil')}
    candidate = c['solution']
    if candidate.has(sp.Integral) or candidate.free_symbols - {x}:
        return {'u': None, 'status': 'not_found', 'checks': [], 'reason': 'kapalı biçimde değil'}
    verdict = verify_exact(r, candidate)
    return {'u': candidate if verdict['status'].startswith('verified') else None, 'candidate': candidate, **verdict}



# ---------- 4a. ÖĞRETİCİ TÜRETMELER (kâğıt; hesap motorunu kullanmaz, yalnız onunla karşılaştırılır) ----------
# Neden: kâğıt sonuçların NEREDEN geldiğini göstermeli. Aşağıdaki yardımcılar her derece ve her problem türü için aynı
# algoritmayı uygular; ürettikleri her ifade hesap motorunun kullandığıyla karşılaştırılır, uyuşmazsa kâğıt üretilmez
# (gösterilen türetme ile yapılan hesap ayrışamaz).
GLOBAL_SHOW = 2       # eleman başına ayrıntı (global N, u_h, R) en fazla bu kadar elemanda; diğerleri aynı yapıda
ROW_SHOW = 9          # test fonksiyonu başına denklem bu kadar satıra kadar tam yazılır


def hat(k):
    return sp.Symbol(f'\\hat{{a}}^{{(e)}}_{{{k}}}')      # e elemanının yerel katsayısı (global a_j ile karıştırılmasın)


def node_value(j):
    return sp.Symbol(f'a_{{{j}}}')


# Şekil fonksiyonları polinom interpolasyonundan: u^e(s) = Σ c_m s^m, u^e(s_k) = â_k → katsayılar → âₖ'ya göre gruplama.
def shape_function_derivation(degree, h, N_engine):
    cs = [sp.Symbol(f'c_{{{m}}}') for m in range(degree + 1)]
    hats = [hat(k + 1) for k in range(degree + 1)]
    nodes = [sp.Rational(k, degree) * h for k in range(degree + 1)]
    poly = sum(c * s ** m for m, c in enumerate(cs))
    blocks = []
    B = blocks.append
    B(f'u^{{e}}(s) = {sp.latex(poly)},\\qquad 0 \\le s \\le h = {sp.latex(h)}\\qquad '
      f'\\text{{({degree}. derece polinom: {degree + 1} katsayı, {degree + 1} düğüm)}}')
    B('s_k = (k - 1)\\,\\frac{h}{p}:\\qquad ' + ',\\ '.join(f's_{{{k + 1}}} = {sp.latex(v)}' for k, v in enumerate(nodes)))
    conds = [sp.expand(poly.subs(s, v)) for v in nodes]
    B('\\begin{aligned} ' + ' \\\\ '.join(f'u^{{e}}(s_{{{k + 1}}}) &= {sp.latex(cond)} = {sp.latex(hats[k])}'
                                          for k, cond in enumerate(conds)) + ' \\end{aligned}\\qquad \\text{(interpolasyon koşulları)}')
    V = sp.Matrix([[v ** m for m in range(degree + 1)] for v in nodes])
    B(f'{tex_matrix(V.tolist())}\\,{tex_column(cs)} = {tex_column(hats)}\\qquad \\mathbf{{V}}\\mathbf{{c}} = \\hat{{\\mathbf{{a}}}}^{{(e)}}\\quad \\text{{(interpolasyon matrisi }} \\mathbf{{V}}\\text{{, }} V_{{km}} = s_k^{{\\,m}}\\text{{)}}')
    sol = sp.solve([sp.Eq(cond, hv) for cond, hv in zip(conds, hats)], cs, dict=True)[0]
    if degree >= 1:
        B(f's_1 = 0 \;\\Rightarrow\; c_0 = {sp.latex(sol[cs[0]])}')
        W = V[1:, 1:]
        rhs = [hv - hats[0] for hv in hats[1:]]
        B(f'{tex_matrix(W.tolist())}\\,{tex_column(cs[1:])} = {tex_column(rhs)}\\qquad \\text{{(kalan {degree} denklem, }} c_0 \\text{{ yerine konmuş)}}')
    B('\\begin{aligned} ' + ' \\\\ '.join(f'{sp.latex(c)} &= {sp.latex(sp.collect(sp.expand(sol[c]), hats))}' for c in cs)
      + ' \\end{aligned}\\qquad \\text{(çözüm)}')
    substituted = sp.expand(poly.subs(sol))
    shapes = [sp.expand(substituted.coeff(hv)) for hv in hats]
    B(f'u^{{e}}(s) = ' + ' + '.join(f'{sp.latex(hv)}\\left({sp.latex(Nk)}\\right)' for hv, Nk in zip(hats, shapes))
      + f' = \\sum_{{k=1}}^{{{degree + 1}}} \\hat{{a}}^{{(e)}}_k\\,N_k(s)\\qquad \\text{{(}}\\hat{{a}}^{{(e)}}_k\\text{{ katsayılarına göre gruplandı)}}')
    for k, (Nk, Ne) in enumerate(zip(shapes, N_engine)):
        if sp.simplify(Nk - Ne) != 0:
            raise RuntimeError(f'Türetilen N_{k + 1} hesapta kullanılanla aynı değil')    # gösterim ile hesap ayrışamaz
    B('\\begin{aligned} ' + ' \\\\ '.join(
        f"N_{{{k + 1}}}(s) &= {sp.latex(Nk)}, &\\qquad N_{{{k + 1}}}'(s) &= {sp.latex(sp.diff(Nk, s))}" for k, Nk in enumerate(shapes))
        + ' \\end{aligned}')
    delta = sp.Matrix(degree + 1, degree + 1, lambda i, j: shapes[i].subs(s, nodes[j]))
    B(f'\\big[N_i(s_j)\\big] = {tex_matrix(delta.tolist())} = \\mathbf{{I}}\\quad (N_i(s_j) = \\delta_{{ij}}),\\qquad '
      f'\\sum_k N_k(s) = {sp.latex(sp.simplify(sum(shapes)))}')
    return blocks


# Elemanın global koordinattaki yaklaşımı: x = x_e + s, â_k = a_{(e−1)p+k}.
def element_approximation(r, e):
    degree, h, x0 = r['degree'], r['h'], r['x0']
    xe = x0 + e * h
    Ns = [sp.expand(Nk.subs(s, x - xe)) for Nk in r['N']]
    nodes = [node_value(e * degree + k + 1) for k in range(degree + 1)]
    return xe, Ns, nodes, sum(a * Nk for a, Nk in zip(nodes, Ns))


# Elemandaki kalıntı: R^{(e)} = L[u_h^{(e)}] − f. Eleman İÇİNDE klasik türevle; 2. mertebede (αu_h')' terimi de dahil.
def element_residual(r, uh):
    if r['order'] == 1:
        R = r['beta'] * sp.diff(uh, x) + r['gamma'] * uh - r['f']
    else:
        R = sp.diff(r['alpha'] * sp.diff(uh, x), x) + r['beta'] * sp.diff(uh, x) + r['gamma'] * uh - r['f']
    syms = sorted(uh.free_symbols - {x}, key=str)
    return sp.collect(sp.expand(R), syms)

# ---------- 4. KÂĞIT ----------
# Sürekli sınıflandırma ve katsayı dejenereliği kâğıtta kanıt düzeyiyle yazılır (sembolik / sayısal / yapılamadı).
def classification_lines(r):
    lines, c, g = [], r.get('classification'), r.get('degeneracy')
    lead = '\\alpha' if r['order'] == 2 else '\\beta'
    if g and g['kind'] != 'uniform':
        where = ',\\ '.join(sp.latex(pt) if g['method'] == 'sembolik' else dec(pt, 6) for pt in g['points'])
        name = 'işaret değiştiren' if g['kind'] == 'sign_changing' else 'dejenere'
        lines.append(f'{lead}(x^*) = 0,\\ x^* \\in \\{{{where}\\}}\\qquad \\text{{({name} katsayı; {g["method"]} tespit; düzgün eliptiklik yok, standart yakınsama kuramı garanti değil)}}')
    if c:
        if c['kind'] == 'unique':
            extra = (';\\ ' + '; '.join(c['assumptions'])) if c['assumptions'] else ''
            lines.append(f'\\text{{Sürekli problem: tek çözüm ({c["method"]}; {c["evidence"]}{extra})}}')
        elif c['kind'] in ('no_solution', 'non_unique'):          # yalnız sayısal kanıtlı olanlar buraya ulaşır
            name = 'çözüm yok' if c['kind'] == 'no_solution' else 'tek değil'
            lines.append(f'\\text{{Sürekli problem: {name} ({c["method"]}; kanıt değil; {c["evidence"]})}}')
        else:
            lines.append(f'\\text{{Sürekli problem sınıflandırılamadı: {c["reason"]}}}')
    return lines


def build_sheet(r):
    """Kâğıt: ders kitabı ve ödevdeki matematiksel sırayla. Güçlü form → gösterim → ağ → şekil fonksiyonu türetmesi →
    yaklaşık çözüm → Galerkin ağırlıklı kalıntı → (yalnız 2. mertebede) zayıf form → Ka = F'nin türetilmesi → eleman
    integralleri → montaj eşlemesi → sınır koşulları ve indirgenmiş sistem → çözüm → yeniden kurulum → karşılaştırma.
    Her ifade hesap motorunun gerçekten yaptığı işlemle aynıdır (türetmeler hesapla otomatik karşılaştırılır)."""
    alpha, beta, gamma, f = r['alpha'], r['beta'], r['gamma'], r['f']
    x0, x1, n, degree, h, dofs, order = r['x0'], r['x1'], r['n'], r['degree'], r['h'], r['dofs'], r['order']
    blocks = []
    add = blocks.append
    D = lambda body: add(f'$${body}$$')
    X0, X1 = sp.latex(x0), sp.latex(x1)
    number = [0]

    def section(title):                       # bölüm numarası probleme göre (zayıf form bölümü yalnız gerektiğinde var)
        number[0] += 1
        add(f'## {number[0]}. {title}')

    def condition_tex(side, point):
        kind, value = r['conditions'][side]
        P = sp.latex(point)
        if kind == 'dirichlet':
            return f'u({P}) = {sp.latex(value)}\\ \\text{{(temel koşul)}}'
        if kind == 'neumann':
            return f"u'({P}) = {sp.latex(value)}\\ \\text{{(doğal koşul: akı }}\\alpha u'\\text{{ verilir)}}"
        if kind == 'both':
            return f"u({P}) = {sp.latex(value)},\\ u'({P}) = {sp.latex(r['slopes'][side])}\\ \\text{{(başlangıç koşulları)}}"
        if r['ivp'] or order != 2:
            return None                             # IVP'nin serbest ucu ve 1. mertebede koşulsuz uç: koşul yok
        if alpha.subs(x, point) == 0:
            return f"\\alpha({P}) = 0\\ \\text{{(verilmedi: akı terimi kendiliğinden sıfır, bağımsız koşul değil)}}"
        return f"\\alpha u'\\big|_{{x={P}}} = 0\\ \\text{{(verilmedi: doğal koşul)}}"

    lhs_tex = tex_equation(alpha, beta, gamma, f).split(' = ')[0]
    add('## Verilenler')
    D(tex_equation(alpha, beta, gamma, f) + f', \\qquad x \\in \\Omega = ({X0},\\, {X1})')
    D(f'\\mathcal{{L}}(u) = {lhs_tex}\\qquad \\text{{({order}. mertebe diferansiyel operatör; güçlü form: }} \\mathcal{{L}}(u) = f\\text{{)}}')
    D(',\\qquad '.join(c for c in (condition_tex('left', x0), condition_tex('right', x1)) if c))
    D(f'n = {n},\\qquad h = \\frac{{x_1 - x_0}}{{n}} = {sp.latex(h)},\\qquad p = {degree}\\ \\text{{({DEGREE_NAMES[degree]} eleman: eleman başına {degree + 1} düğüm)}},'
      f'\\qquad N = np + 1 = {dofs}\\ \\text{{(toplam düğüm = bilinmeyen sayısı)}}')
    method_text = ('Bubnov–Galerkin' if not r['ivp'] else 'Petrov–Galerkin (test fonksiyonları akısı bilinmeyen uçta sıfır)')
    D(f'\\text{{Yöntem: {method_text}}}')
    for line in classification_lines(r):
        D(line)

    section('Gösterim')
    D('\\begin{aligned} u(x)&:\\ \\text{tam çözüm (biliniyorsa)} & u_h(x)&:\\ \\text{sonlu eleman yaklaşımı} \\\\ '
      '\\phi_j(x)&:\\ \\text{global yaklaşım (şekil) fonksiyonu} & N_k(s),\\ N_k^{(e)}(x)&:\\ \\text{yerel şekil fonksiyonu, } e \\text{ elemanında} \\\\ '
      'a_j&:\\ \\text{bilinmeyen global düğüm katsayısı} & \\hat{a}^{(e)}_k&:\\ e \\text{ elemanının yerel katsayısı} \\\\ '
      'c_m&:\\ \\text{eleman polinomunun katsayısı} & R(u_h, x)&:\\ \\text{kalıntı} \\\\ '
      'w_i(x)&:\\ \\text{ağırlık (test) fonksiyonu} & K_{ij},\\ F_i&:\\ \\text{sistem matrisi ve yük vektörü girdileri} \\end{aligned}')
    D('\\phi_j(x) = \\begin{cases} N_k^{(e)}(x), & x \\in \\Omega_e,\\ j = (e-1)p + k \\\\ 0, & \\text{diğer elemanlarda} \\end{cases}'
      '\\qquad \\text{(ders kitabındaki } \\phi_j \\text{, sonlu eleman interpolasyonunda parça parça } N_k^{(e)} \\text{ olur)}')

    section('Ağ')
    nodes = [x0 + sp.Rational(i, degree) * h for i in range(dofs)]
    D('\\big[x_j\\big] = ' + tex_row(nodes) + f'\\qquad \\Omega_e = [x_e,\\ x_e + h],\\ e = 1, \\dots, {n}')
    shown = list(range(n)) if n <= 6 else [0, 1, None, n - 1]
    rows = ' \\\\ '.join('\\vdots & \\vdots & \\vdots' if e is None else
                         f'{e + 1} & [{sp.latex(x0 + e * h)},\\ {sp.latex(x0 + (e + 1) * h)}] & '
                         + ',\\ '.join(f'{e * degree + k + 1}' for k in range(degree + 1)) for e in shown)
    D(f'\\begin{{array}}{{c|c|c}} e & \\Omega_e & \\text{{yerel }} k = 1, \\dots, {degree + 1} \;\\mapsto\; \\text{{global }} j = (e-1)p + k \\\\ \\hline {rows} \\end{{array}}')

    section('Şekil fonksiyonları (polinom interpolasyonundan)')
    D(f's = x - x_e \\in [0, h]\\qquad \\text{{(yerel koordinat: her elemanda aynı türetme)}}')
    for block in shape_function_derivation(degree, h, r['N']):
        D(block)
    show = r['detailed'][:GLOBAL_SHOW]
    D('x = x_e + s,\\quad dx = ds,\\qquad N_k^{(e)}(x) = N_k(x - x_e),\\qquad \\hat{a}^{(e)}_k = a_{(e-1)p+k}\\qquad \\text{(yerel → global)}')
    for e in show:
        xe, Ns, nodes_e, _ = element_approximation(r, e)
        D(f'e = {e + 1}:\\qquad ' + ',\\quad '.join(f'N_{{{k + 1}}}^{{({e + 1})}}(x) = {sp.latex(Nk)}' for k, Nk in enumerate(Ns))
          + ',\\qquad ' + ',\\ '.join(f'\\hat{{a}}^{{({e + 1})}}_{{{k + 1}}} = {sp.latex(a)}' for k, a in enumerate(nodes_e)))

    section('Yaklaşık çözüm')
    D(f'u_h(x) = \\sum_{{j=1}}^{{{dofs}}} a_j\\,\\phi_j(x)\\qquad \\text{{(}} a_j \\text{{: belirlenecek bilinmeyenler)}}')
    D('\\phi_j(x_i) = \\delta_{ij} \;\\Rightarrow\; u_h(x_i) = a_i\\qquad \\text{(katsayılar düğüm değerleridir)}')
    D('x \\in \\Omega_e:\\quad u_h(x) = u_h^{(e)}(x) = \\sum_{k=1}^{p+1} \\hat{a}^{(e)}_k\\,N_k^{(e)}(x)')
    for e in show:
        xe, Ns, nodes_e, uh = element_approximation(r, e)
        D(f'u_h^{{({e + 1})}}(x) = ' + ' + '.join(f'{sp.latex(a)}\\left({sp.latex(Nk)}\\right)' for a, Nk in zip(nodes_e, Ns))
          + f',\\qquad x \\in [{sp.latex(xe)}, {sp.latex(xe + h)}]')
    if n > len(show):
        D('\\text{diğer elemanlar aynı yapıda}')

    # --- Galerkin ağırlıklı kalıntı
    section('Galerkin Weighted-Residual Formulation (Galerkin ağırlıklı kalıntı)')
    D(f'R(u_h, x) = \\mathcal{{L}}(u_h) - f = \\left[{lhs_tex}\\right]_{{u = u_h}}' + ('' if f == 0 else f' - \\left({sp.latex(f)}\\right)'))
    residuals = {}
    for e in show:
        xe, Ns, nodes_e, uh = element_approximation(r, e)
        residuals[e] = element_residual(r, uh)
        D(f'R(u_h, x)\\big|_{{\\Omega_{{{e + 1}}}}} = {sp.latex(residuals[e])}')
    D(f'\\int_{{{X0}}}^{{{X1}}} w_i(x)\\,R(u_h, x)\\,dx = 0\\qquad \\text{{(kalıntı, ağırlık fonksiyonlarına göre ortalamada sıfır)}}')
    D('w_i(x) = \\phi_i(x)\\qquad \\text{(' + ('Bubnov–Galerkin: ağırlık fonksiyonları yaklaşım fonksiyonlarıyla aynı' if not r['ivp']
      else 'Petrov–Galerkin: ağırlıklar yaklaşım fonksiyonlarından seçilir, test uzayı farklı') + ')}')
    tests_used = list(r['rows'])
    reasons = []
    for i in range(dofs):
        if i in tests_used:
            continue
        P = 'x_0' if i == 0 else 'x_1'
        if r['ivp']:
            reasons.append(f'w_{{{i + 1}}} = \\phi_{{{i + 1}}}\\ \\text{{çıkarıldı: }} {P}\\text{{ ucunda akı bilinmiyor, test fonksiyonu orada sıfır seçilir}}')
        else:
            reasons.append(f'w_{{{i + 1}}} = \\phi_{{{i + 1}}}\\ \\text{{çıkarıldı: }} u({P})\\ \\text{{verildi (temel koşul) }}\\Rightarrow a_{{{i + 1}}}\\ \\text{{bilinen; test uzayı: }} w({P}) = 0')
    D(f'w_i \\in \\{{\\phi_i\\}},\\quad {index_range(tests_used[0] + 1, tests_used[-1] + 1)}\\qquad \\text{{({len(tests_used)} denklem, {len(r["free"])} bilinmeyen)}}')
    for reason in reasons:
        D(reason)
    a_syms = [node_value(j + 1) for j in range(dofs)]

    def matrix_definitions():
        D(f'\\mathbf{{a}} = \\begin{{bmatrix}} a_1 & \\cdots & a_{{{dofs}}} \\end{{bmatrix}}^{{\\mathsf{{T}}}} \\in \\mathbb{{R}}^{{{dofs}}},\\quad '
          f'\\mathbf{{K}} \\in \\mathbb{{R}}^{{{dofs}\\times{dofs}}},\\quad \\mathbf{{F}} \\in \\mathbb{{R}}^{{{dofs}}}\\qquad '
          '\\text{(satır } i \\leftrightarrow w_i = \\phi_i\\text{, sütun } j \\leftrightarrow a_j\\text{; } K_{ij}\\text{: } a_j\\text{\'nin katsayısı)}')
        D('\\sum_{j} K_{ij}\\,a_j = F_i \;\\Longleftrightarrow\; \\mathbf{K}\\,\\mathbf{a} = \\mathbf{F}')

    if order == 1:
        D('\\text{Birinci mertebe: } R \\text{ yalnız } u_h\' \\text{ içerir, } C^0 \\text{ yaklaşımda integrallenebilir} \\Rightarrow \\text{kısmi integrasyon gerekmez}')
        if dofs <= ROW_SHOW:
            for i in tests_used:
                parts, total = [], 0
                for e in range(n):
                    kloc = i - e * degree
                    if not 0 <= kloc <= degree:
                        continue
                    xe, Ns, nodes_e, uh = element_approximation(r, e)
                    Re = residuals[e] if e in residuals else element_residual(r, uh)
                    parts.append(f'\\int_{{{sp.latex(xe)}}}^{{{sp.latex(xe + h)}}}\\left({sp.latex(Ns[kloc])}\\right)\\left[{sp.latex(Re)}\\right]dx')
                    total += sp.integrate(sp.expand(Ns[kloc] * Re), (x, xe, xe + h))
                total = sp.collect(sp.expand(total), a_syms)
                row = sp.expand(sum(r['K'][i, j] * a_syms[j] for j in range(dofs)) - r['F_l'][i])
                if total.has(sp.Integral) or (not r['K'].has(sp.Float) and sp.simplify(total - row) != 0):
                    raise RuntimeError(f'Galerkin denklemi i={i + 1}, K satırıyla tutmuyor')        # gösterim = hesap
                D(f'w_{{{i + 1}}}:\\quad ' + ' + '.join(parts) + f' = {sp.latex(total)} = 0')
        k_int = tex_terms([(beta, "\\phi_i\\,\\phi_j'"), (gamma, '\\phi_i\\,\\phi_j')])
        f_int = '0' if f == 0 else f'\\int_{{{X0}}}^{{{X1}}} \\phi_i\\left({sp.latex(f)}\\right)dx'
        D(f'u_h = \\sum_j a_j\\phi_j:\\quad \\int_{{{X0}}}^{{{X1}}} \\phi_i\\,R\\,dx = \\sum_j a_j \\underbrace{{\\int_{{{X0}}}^{{{X1}}} \\left({k_int}\\right)dx}}_{{K_{{ij}}}}'
          f' - \\underbrace{{{f_int}}}_{{F_i}} = 0')
        matrix_definitions()
    else:
        D("\\text{İkinci mertebe: } R \\ni (\\alpha u_h')' \\text{; } C^0 \\text{ yaklaşımda } u_h' \\text{ düğümlerde süreksiz, } (\\alpha u_h')' \\text{ orada klasik anlamda tanımsız}"
          "\\Rightarrow \\text{türev ağırlık fonksiyonuna aktarılır (bir sonraki bölüm)}")
        D(f"\\int_{{{X0}}}^{{{X1}}} w_i\\left[{lhs_tex.replace('u', 'u_h')}" + ('' if f == 0 else f' - \\left({sp.latex(f)}\\right)') + '\\right]dx = 0')

        section('Weak Formulation (Integration by Parts)')
        D(f"\\int_{{{X0}}}^{{{X1}}} w_i\\,(\\alpha u_h')'\\,dx = \\Big[w_i\\,\\alpha u_h'\\Big]_{{{X0}}}^{{{X1}}} - \\int_{{{X0}}}^{{{X1}}} w_i'\\,\\alpha u_h'\\,dx"
          "\\qquad \\text{(kısmi integrasyon)}")
        D("\\text{Sonuç: } u_h\\text{'den yalnız birinci türev istenir (türev mertebesi 2 → 1); } C^0 \\text{ yaklaşım yeterli; uçlarda akı } Q = \\alpha u' \\text{ ortaya çıkar}")
        weak_left = tex_terms([(alpha, "w_i'u_h'"), (-beta, "w_i\\,u_h'"), (-gamma, 'w_i\\,u_h')])
        load = '' if f == 0 else f" - \\int_{{{X0}}}^{{{X1}}} w_i\\left({sp.latex(f)}\\right)dx"
        D(f"\\int_{{{X0}}}^{{{X1}}}\\left({weak_left}\\right)dx = w_i({X1})\\,Q_1 - w_i({X0})\\,Q_0{load},\\qquad Q = \\alpha u'")
        for side, point, name in (('left', x0, 'Q_0'), ('right', x1, 'Q_1')):
            kind, value = r['conditions'][side]
            P = sp.latex(point)
            if kind == 'dirichlet':
                why = f'u({P})\\ \\text{{verildi: test uzayında }} w_i({P}) = 0 \\Rightarrow \\text{{terim kullanılan denklemlerde yok; }} {name} \\text{{ bilinmeyen, çözümden sonra hesaplanır}}'
            elif kind == 'neumann':
                why = f"{name} = \\alpha({P})\\,u'({P}) = {tex_value(r['Q'][side])}\\ \\text{{(doğal koşul: bilinen akı, yük vektörüne girer)}}"
            elif kind == 'both':
                why = f"{name} = \\alpha({P})\\,u'({P}) = {tex_value(r['Q'][side])}\\ \\text{{(başlangıç koşulu: akı bilinir)}}"
            elif r['ivp']:
                why = f'{name}\\ \\text{{bilinmiyor: test uzayında }} w_i({P}) = 0 \\text{{ seçilir (Petrov–Galerkin), terim düşer}}'
            elif alpha.subs(x, point) == 0:
                why = f'\\alpha({P}) = 0 \\Rightarrow {name} = 0\\ \\text{{(dejenere uç)}}'
            else:
                why = f'{name} = 0\\ \\text{{(koşul verilmedi: doğal koşul }} \\alpha u\' = 0\\text{{)}}'
            D(f'x = {P}:\\quad ' + why)
        k_int = tex_terms([(alpha, "\\phi_i'\\phi_j'"), (-beta, "\\phi_i\\,\\phi_j'"), (-gamma, '\\phi_i\\,\\phi_j')])
        f_int = '' if f == 0 else f' - \\int_{{{X0}}}^{{{X1}}} \\phi_i\\left({sp.latex(f)}\\right)dx'
        D(f'u_h = \\sum_j a_j\\phi_j,\\ w_i = \\phi_i:\\quad \\sum_j a_j \\underbrace{{\\int_{{{X0}}}^{{{X1}}} \\left({k_int}\\right)dx}}_{{K_{{ij}}}}'
          f' = \\underbrace{{\\phi_i({X1})\\,Q_1 - \\phi_i({X0})\\,Q_0{f_int}}}_{{F_i = F_{{b,i}} + F_{{l,i}}}}')
        matrix_definitions()

    # --- eleman integralleri
    section('Eleman matrisleri ve yük vektörleri')
    D(f'K_{{ij}} = \\sum_e K^{{(e)}}_{{ij}},\\qquad K^{{(e)}}_{{ij}} = \\int_{{\\Omega_e}} (\\cdots)\\,dx \\text{{ (yukarıdaki integrand, }} \\phi \\to N^{{(e)}}\\text{{)}},\\qquad '
      f'i, j = 1, \\dots, {degree + 1}\\ \\text{{(yerel)}},\\qquad \\mathbf{{K}}^{{(e)}} \\in \\mathbb{{R}}^{{{degree + 1}\\times{degree + 1}}},\\ \\mathbf{{F}}^{{(e)}} \\in \\mathbb{{R}}^{{{degree + 1}}}')
    integrand_K = (tex_terms([(alpha, "N_i'N_j'"), (-beta, "N_i\\,N_j'"), (-gamma, 'N_i\\,N_j')]) if order == 2
                   else tex_terms([(beta, "N_i\\,N_j'"), (gamma, 'N_i\\,N_j')]))
    D(f'K^{{(e)}}_{{ij}} = \\int_0^h\\left({integrand_K}\\right)ds,\\qquad F^{{(e)}}_{{l,i}} = {"-" if order == 2 else ""}\\int_0^h N_i\\,f(x_e + s)\\,ds,\\qquad x = x_e + s')
    if any('\\approx' in line for record in r['elements'] for line in record['K_lines'] + record['F_lines']):
        D('\\approx:\\ \\text{kapalı biçim bulunamadı; sayısal integral (mpmath.quad, uyarlamalı tanh–sinh kuralı, 15 basamak)}')
    else:
        D('\\text{İntegraller SymPy ile tam (sembolik) alındı.}')
    if r['uniform_K'] and n > 1:
        D('\\alpha, \\beta, \\gamma \\text{ sabit ve } h \\text{ eşit} \\Rightarrow x = x_e + s \\text{ ile integrand } e\\text{\'den bağımsız} \\Rightarrow \\mathbf{K}^{(e)} = \\mathbf{K}^{(1)}')
    if r['uniform_F'] and n > 1:
        D('f \\text{ sabit} \\Rightarrow \\mathbf{F}^{(e)}_l = \\mathbf{F}^{(1)}_l')
    for record in r['elements']:
        e = record['e']
        if e not in r['detailed']:
            continue
        has_K, has_F = not (r['uniform_K'] and e > 0), not (r['uniform_F'] and e > 0)
        if not (has_K or has_F):
            continue
        D(f'e = {e + 1}:\\qquad x_e = {sp.latex(record["xe"])},\\qquad x = {sp.latex(record["xe"] + s)}')
        if has_K:
            for line in record['K_lines']:
                D(line)
            D(f'\\mathbf{{K}}^{{({e + 1})}} = {tex_matrix(record["K"].tolist())}')
        if has_F:
            for line in record['F_lines']:
                D(line)
            D(f'\\mathbf{{F}}_l^{{({e + 1})}} = {tex_column(list(record["F"]))}')
    if n > DETAIL_LIMIT and not (r['uniform_K'] and r['uniform_F']):
        D(f'e = 3, \\dots, {n - 1}:\\quad \\text{{aynı integraller; değerleri montajda}}')
    if order == 2 and not alpha.free_symbols and alpha < 0:
        first_K, first_F = r['elements'][0]['K'], r['elements'][0]['F']
        D(f'\\alpha = {sp.latex(alpha)} < 0:\\qquad \\tilde{{\\mathbf{{K}}}}^{{(1)}} = -\\mathbf{{K}}^{{(1)}} = {tex_matrix((-first_K).tolist())},'
          f'\\qquad \\tilde{{\\mathbf{{F}}}}_l^{{(1)}} = -\\mathbf{{F}}_l^{{(1)}} = {tex_column(list(-first_F))}')
        D('\\tilde{\\mathbf{K}}\\,\\mathbf{a} = \\tilde{\\mathbf{F}} \\iff (-\\mathbf{K})\\,\\mathbf{a} = -\\mathbf{F} \\iff '
          '\\mathbf{K}\\,\\mathbf{a} = \\mathbf{F}\\qquad \\text{(aynı }\\mathbf{a}\\text{; standart pozitif rijitlik işaretiyle)}')

    # --- montaj
    section('Montaj')
    D('K_{IJ} = \\sum_{e} K^{(e)}_{ij},\\qquad F_I = \\sum_e F^{(e)}_i,\\qquad I = (e-1)p + i,\\ J = (e-1)p + j\\qquad \\text{(yerel } i,j \\to \\text{ global } I,J\\text{)}')
    if dofs <= ROW_SHOW:
        contrib = {}
        for record in r['elements']:
            base = record['e'] * degree
            for i in range(degree + 1):
                for j in range(degree + 1):
                    contrib.setdefault((base + i, base + j), []).append((record['e'], i, j))
        shared = [(IJ, c) for IJ, c in sorted(contrib.items()) if len(c) > 1]
        for (I, J), c in shared:
            terms = ' + '.join(f'K^{{({e + 1})}}_{{{i + 1}{j + 1}}}' for e, i, j in c)
            values = ' + '.join(f'\\left({tex(r["elements"][e]["K"][i, j])}\\right)' for e, i, j in c)
            D(f'K_{{{I + 1}{J + 1}}} = {terms} = {values} = {tex(r["K"][I, J])}\\qquad \\text{{(ortak düğüm {I + 1})}}' if I == J else
              f'K_{{{I + 1}{J + 1}}} = {terms} = {values} = {tex(r["K"][I, J])}')
        fshared = {}
        for record in r['elements']:
            for i in range(degree + 1):
                fshared.setdefault(record['e'] * degree + i, []).append((record['e'], i))
        for I, c in sorted(fshared.items()):
            if len(c) > 1:
                D(f'F_{{l,{I + 1}}} = ' + ' + '.join(f'F^{{({e + 1})}}_{{l,{i + 1}}}' for e, i in c) + ' = '
                  + ' + '.join(f'\\left({tex(r["elements"][e]["F"][i])}\\right)' for e, i in c) + f' = {tex(r["F_l"][I])}')
        D('\\text{diğer girdiler tek bir elemandan gelir: } K_{IJ} = K^{(e)}_{ij}' + ',\\ \\text{boş girdiler: ortak elemanı olmayan düğümler}')
    D('\\mathbf{K} = ' + tex_matrix(r['K'].tolist()) + ',\\qquad \\mathbf{F}_l = ' + tex_column(list(r['F_l'])))
    if order == 2:
        D("\\mathbf{F}_b = \\begin{bmatrix} -Q_0 & 0 & \\cdots & 0 & Q_1 \\end{bmatrix}^{\\mathsf{T}}\\qquad \\text{(yalnız uç düğümler; iç düğümlerde } \\phi_i(x_0) = \\phi_i(x_1) = 0\\text{)}")
        if dofs <= ROW_SHOW:
            Qs = {0: sp.Symbol('Q_0'), dofs - 1: sp.Symbol('Q_1')}
            for i in r['rows']:
                lhs = sp.expand(sum(r['K'][i, j] * a_syms[j] for j in range(dofs)))
                Fb = r['F_b'][i] if (i not in Qs or r['Q']['left' if i == 0 else 'right'] is not None) else (-Qs[i] if i == 0 else Qs[i])
                D(f'w_{{{i + 1}}}:\\quad {sp.latex(lhs)} = {sp.latex(sp.expand(r["F_l"][i] + Fb))}')

    # --- sınır koşulları
    section('Sınır koşulları ve indirgenmiş sistem')
    known = []
    for side, index, name in (('left', 0, 'Q_0'), ('right', dofs - 1, 'Q_1')):
        if index in r['fixed']:
            known.append(f'a_{{{index + 1}}} = {sp.latex(r["fixed"][index])}\\ \\text{{(temel)}}')
        if order == 2 and r['Q'][side] is not None:
            known.append(f'{name} = {tex_value(r["Q"][side])}\\ \\text{{(doğal)}}')
    if known:
        D(',\\qquad '.join(known))
    removed_rows = [i + 1 for i in range(dofs) if i not in r['rows']]
    removed_cols = [i + 1 for i in r['fixed']]
    if removed_rows or removed_cols:
        D(('\\text{çıkarılan satırlar: } ' + ', '.join(str(i) for i in removed_rows) + '\\ (w_i \\text{ test uzayında değil})' if removed_rows else '')
          + ('\\qquad ' if removed_rows and removed_cols else '')
          + ('\\text{çıkarılan sütunlar: } ' + ', '.join(str(j) for j in removed_cols) + '\\ (a_j \\text{ bilinen} \\Rightarrow K_{ij}a_j \\text{ sağ tarafa: } F_{\\mathrm{red},i} = F_i - \\sum_{j\\ \\text{bilinen}} K_{ij}a_j)' if removed_cols else ''))
    if r['free']:
        unknowns = [sp.Symbol(f'a_{{{i + 1}}}') for i in r['free']]
        D(f'\\mathbf{{K}}_{{\\mathrm{{red}}}}\\,\\mathbf{{a}}_{{\\mathrm{{bilinmeyen}}}} = \\mathbf{{F}}_{{\\mathrm{{red}}}}:\\qquad {tex_matrix(r["K_ff"].tolist())}\\,{tex_column(unknowns)} = {tex_column(list(r["F_f"]))}')

    # --- çözüm
    section('Çözüm')
    D('\\mathbf{a} = ' + tex_row(r['a']))
    if r['exact'] and any(not sp.sympify(v).is_Rational for v in r['a']):
        D('\\mathbf{a} \\approx ' + tex_row(r['a'], lambda v: dec(sp.sympify(v).evalf())))
    for side, name, point in (('left', 'Q_0', x0), ('right', 'Q_1', x1)):
        if side in r['reactions']:
            flux, a_end = r['reactions'][side], alpha.subs(x, point)
            slope = f",\\qquad u_h'({sp.latex(point)}) = \\frac{{{name}}}{{\\alpha({sp.latex(point)})}} = {tex_value(sp.simplify(flux / a_end))}" if a_end != 0 else ''
            D(f'{name} = {tex_value(flux)}{slope}\\qquad \\text{{(çıkarılan satırın denkleminden)}}')
    if r['first_order']:
        fo = r['first_order']
        rd, d = fo['residual'], fo['row']
        D(f"\\text{{çıkarılan denklem }} w_{{{d}}}\\text{{: kalıntısı }} r_{{{d}}} = \\sum_j K_{{{d}j}}\\,a_j - F_{{{d}}} = {tex(rd) if fo['exact'] else dec(rd)}\\qquad \\text{{(otomatik sağlanmaz)}}")
        yb_tex = tex(fo['yb']) if fo['exact'] and fo['yb'] is not None else dec(fo['yb'])
        scaled = '' if fo['exact'] else f",\\quad \\frac{{|\\mathbf{{y}}^{{\\mathsf{{T}}}}\\mathbf{{b}}|}}{{\\|\\mathbf{{y}}\\|\\,\\|\\mathbf{{b}}\\|}} = {dec(fo['scaled'])}"
        D(f"\\mathbf{{y}}^{{\\mathsf{{T}}}}\\mathbf{{A}} = 0,\\ \\mathbf{{y}}^{{\\mathsf{{T}}}}\\mathbf{{b}} = {yb_tex}{scaled}"
          f"\\qquad (\\mathbf{{A}}\\mathbf{{a}} = \\mathbf{{b}}:\\ {dofs}\\times{dofs - 1}\\ \\text{{tam sistem, sol sıfır uzayı ile uyumluluk}})")
        if fo['compatible']:
            D('\\text{tam sistem uyumlu: çıkarılan denklem de sağlanıyor; satır seçimi sonucu değiştirmez}')
        elif fo['compatible'] is False:
            D('\\text{tam sistem uyumsuz: başka bir satırın çıkarılması (başka test uzayı) başka } \\mathbf{a}\\text{ verir}')
    if r['ivp']:
        side = 'left' if r['conditions']['left'][0] == 'both' else 'right'
        e, s_at, point = (0, 0, x0) if side == 'left' else (n - 1, h, x1)
        slope_h = sp.simplify(sum(sp.sympify(r['a'][e * degree + k]) * r['dN'][k].subs(s, s_at) for k in range(degree + 1)))
        D(f"u_h'({sp.latex(point)}) = {tex_value(slope_h)},\\qquad u'({sp.latex(point)}) = {sp.latex(r['slopes'][side])}"
          f"\\qquad \\text{{(u' koşulu akı üzerinden, zayıf biçimde sağlanır)}}")

    # --- yeniden kurulum
    section('Yaklaşık çözümün yeniden kurulması')
    D('u_h^{(e)}(x) = \\sum_{k=1}^{p+1} \\hat{a}^{(e)}_k\\,N_k^{(e)}(x),\\qquad \\hat{a}^{(e)}_k = a_{(e-1)p+k}')
    for e in r['detailed']:
        xe = x0 + e * h
        piece = sp.expand(sum(r['a'][e * degree + k] * r['N'][k] for k in range(degree + 1)).subs(s, x - xe))
        if e in r['detailed'][:GLOBAL_SHOW]:
            _, Ns, _, _ = element_approximation(r, e)
            sub = ' + '.join(f'\\left({tex(sp.sympify(r["a"][e * degree + k]))}\\right)\\left({sp.latex(Nk)}\\right)' for k, Nk in enumerate(Ns))
            D(f'u_h^{{({e + 1})}}(x) = {sub} = {tex(piece)},\\qquad {sp.latex(xe)} \\le x \\le {sp.latex(xe + h)}')
        else:
            D(f'u_h^{{({e + 1})}}(x) = {tex(piece)},\\qquad {sp.latex(xe)} \\le x \\le {sp.latex(xe + h)}')
    if n > DETAIL_LIMIT:
        D('\\text{diğer elemanlarda aynı biçimde}')

    # --- karşılaştırma
    verdict = exact_solution(r)
    solution = verdict['u']
    level = {'verified_symbolic': 'sembolik doğrulandı', 'verified_numeric': 'sayısal doğrulandı (50 basamak; kanıt değil)'}
    if verdict['status'] == 'rejected':
        section('Tam çözümle karşılaştırma')
        failed = next(c for c in verdict['checks'] if not c[1])
        D(f"\\text{{aday }} u(x) = {sp.latex(verdict['candidate'])}\\qquad \\text{{reddedildi: }} {failed[0][1]}\\ \\text{{sağlanmıyor}}\\qquad \\text{{(hata hesaplanmadı)}}")
    if solution is not None:
        section('Tam çözümle karşılaştırma')
        D(f'u(x) = {sp.latex(solution)}\\qquad \\text{{({level[verdict["status"]]}: denklem ve koşullar denetlendi)}}')
        D(',\\quad '.join(f'{name[1]}\\ \\checkmark\\ \\text{{({how})}}' for name, ok, how in verdict['checks']))
        picks = list(range(dofs)) if dofs <= 13 else sorted(set(round(i * (dofs - 1) / 10) for i in range(11)))
        exact_or_dec = lambda v: sp.latex(v) if sp.sympify(v).is_Rational else dec(sp.sympify(v).evalf())
        rows = ' \\\\ '.join(
            f'{i + 1} & {tex(nodes[i])} & {exact_or_dec(r["a"][i])} & {exact_or_dec(sp.nsimplify(solution.subs(x, nodes[i])) if solution.subs(x, nodes[i]).is_Rational else solution.subs(x, nodes[i]))} & '
            f'{exact_or_dec(sp.Abs(sp.sympify(r["a"][i]) - solution.subs(x, nodes[i])))}' for i in picks)
        D(f'\\begin{{array}}{{c|c|c|c|c}} j & x_j & a_j = u_h(x_j) & u(x_j) & |u - u_h| \\\\ \\hline {rows} \\end{{array}}\\qquad \\text{{(düğüm hatası)}}')
        em = error_measures(r, solution)
        mids = ' \\\\ '.join(f'{tex(xm)} & {dec(um)} & {dec(uhm)} & {dec(um - uhm)}' for xm, um, uhm in em['midpoints'])
        D(f'\\begin{{array}}{{c|c|c|c}} x & u(x) & u_h(x) & u - u_h \\\\ \\hline {mids} \\end{{array}}\\qquad \\text{{(noktasal hata, her elemanda ilk iki düğümün ortası }} x_e + h/2p\\text{{)}}')
        D(f'\\max_j |u(x_j) - a_j| = {dec(em["nodal"])}\\ \\text{{(düğüm)}},\\qquad \\max_x |u - u_h| \\approx {dec(em["pointwise"])}\\ \\text{{(noktasal, eleman başına 200 örnek)}}')
        D(f'\\|u - u_h\\|_{{L^2}} = \\left(\\int_{{{X0}}}^{{{X1}}}(u - u_h)^2\\,dx\\right)^{{1/2}} \\approx {dec(em["l2"])},\\qquad '
          f'|u - u_h|_{{H^1}} = \\left(\\int_{{{X0}}}^{{{X1}}}(u\' - u_h\')^2\\,dx\\right)^{{1/2}} \\approx {dec(em["h1"])}')
    r['exact_verdict'] = verdict
    return '\n\n'.join(blocks), solution


def error_measures(r, solution):
    """Kâğıttaki hata ölçüleri, GERÇEK FEM yeniden kurulumuyla (eleman eleman polinom): düğüm hatası, noktasal hata (eleman
    orta noktaları ve eleman başına 200 örnekte en büyük), L² ve H¹ yarı-normu (mpmath.quad her elemanda ayrı; u_h eleman
    sınırında kırılır). Noktasal en büyük değer örneklemeye dayanır (alt tahmin)."""
    degree, h, x0, n = r['degree'], r['h'], r['x0'], r['n']
    exact_f = sp.lambdify(x, solution, 'mpmath')
    exact_d = sp.lambdify(x, sp.diff(solution, x), 'mpmath')
    nodes = [x0 + sp.Rational(i, degree) * h for i in range(r['dofs'])]
    nodal = max(abs(mpmath.mpf(sp.N(sp.sympify(a) - solution.subs(x, xi), 20))) for a, xi in zip(r['a'], nodes))
    l2 = h1 = pointwise = mpmath.mpf(0)
    midpoints = []
    for e in range(n):
        xe = x0 + e * h
        poly = sum(sp.sympify(r['a'][e * degree + k]) * r['N'][k] for k in range(degree + 1))
        piece, dpiece = sp.lambdify(s, poly, 'mpmath'), sp.lambdify(s, sp.diff(poly, s), 'mpmath')
        l2 += mpmath.quad(lambda t: (exact_f(float(xe) + t) - piece(t)) ** 2, [0, float(h)])
        h1 += mpmath.quad(lambda t: (exact_d(float(xe) + t) - dpiece(t)) ** 2, [0, float(h)])
        pointwise = max([pointwise] + [abs(exact_f(float(xe) + float(h) * k / 200) - piece(float(h) * k / 200)) for k in range(201)])
        if e < 6:                                   # ilk iki düğüm arasının ortası: hiçbir p için düğüm değil
            xm = xe + h / (2 * degree)
            midpoints.append((xm, float(exact_f(float(xm))), float(piece(float(h) / (2 * degree)))))
    return {'nodal': float(nodal), 'l2': float(mpmath.sqrt(l2)), 'h1': float(mpmath.sqrt(h1)), 'pointwise': float(pointwise), 'midpoints': midpoints}


# ---------- 5. LATEX BELGESİ VE GİRİŞ ----------
def build_tex(text):
    """Aynı bloklardan derlenebilir belge. \\fit: satıra sığmayan formülü satır genişliğine küçültür (kırpılma olmaz)."""
    body = []
    for block in text.split('\n\n'):
        heading = re.match(r'#+\s+(.*)', block)
        display = re.fullmatch(r'\$\$(.*)\$\$', block, re.S)
        body.append(f'\\section*{{{heading.group(1)}}}' if heading else f'\\[ \\fit{{{display.group(1)}}} \\]' if display else block)
    return '\n'.join([
        '\\documentclass[11pt,a4paper]{article}',
        '\\usepackage[utf8]{inputenc}', '\\usepackage[T1]{fontenc}', '\\usepackage{lmodern,amsmath,amssymb,graphicx}',
        '\\usepackage[margin=18mm]{geometry}', '\\setlength{\\parindent}{0pt}', '\\setcounter{MaxMatrixCols}{20}',
        '\\newcommand{\\fit}[1]{\\resizebox{\\ifdim\\width>\\linewidth\\linewidth\\else\\width\\fi}{!}{$\\displaystyle #1$}}',
        '\\begin{document}', '{\\LARGE\\textbf{Sonlu Eleman Çözüm Kâğıdı}}\\par\\smallskip\\hrule',
        '\n\n'.join(body), '\\end{document}'])


def exact_for_browser(solution, x0, x1):
    """Analitik çözümü tarayıcının ayrıştırıcısının okuyabileceği metne çevirir (u ve u').
    Neden: grafiklerdeki hata, aynı yöntemin ince ağına göre değil, BAĞIMSIZ analitik çözüme göre ölçülebilsin.
    Ondalıklar kesre çevrilir ('1e-5' gibi bilimsel gösterim ayrıştırıcıda e sabiti sanılmasın); tarayıcı metni
    okuduktan sonra 'samples' değerleriyle karşılaştırır, tutmazsa kullanmaz."""
    if solution is None or solution.has(sp.I):
        return None
    solution = solution.xreplace({number: sp.Rational(str(number)) for number in solution.atoms(sp.Float)})
    derivative = sp.diff(solution, x)
    text = lambda expr: sp.sstr(expr).replace('**', '^')
    points = [x0 + (x1 - x0) * sp.Rational(k, 8) for k in range(9)]
    try:
        samples = [[float(pt), float(solution.subs(x, pt).evalf(20)), float(derivative.subs(x, pt).evalf(20))] for pt in points]
    except (TypeError, ValueError):
        return None
    return {'u': text(solution), 'du': text(derivative), 'tex': sp.latex(solution), 'samples': samples}


# Analitik doğrulamanın JSON hali: tarayıcı kanıt düzeyini ve ret nedenini gösterebilsin.
def verdict_payload(verdict):
    return {'exact_status': verdict['status'], 'exact_reason': verdict['reason'],
            'exact_checks': [[name[0], ok, how] for name, ok, how in verdict['checks']],
            'exact_candidate': sp.sstr(verdict['candidate']) if verdict.get('candidate') is not None else None}


# Sınıflandırma ve dejenerelik bilgisinin JSON'a çevrilmiş hali (tarayıcının teşhis kutusu için).
def classification_payload(classification, degeneracy):
    out = {'degeneracy': None, 'classification': None}
    if degeneracy:
        out['degeneracy'] = {'kind': degeneracy['kind'], 'points': [float(pt) for pt in degeneracy['points']], 'method': degeneracy['method']}
    if classification:
        out['classification'] = {k: (sp.sstr(v) if isinstance(v, sp.Basic) else v) for k, v in classification.items()
                                 if k not in ('general', 'solution')}
    return out


def make_sheet(data, symbolic=True):
    result = solve_problem(data, symbolic, classify=True)
    text, solution = build_sheet(result)
    # K: tarayıcının montajladığı matrisle girdi girdi karşılaştırılsın diye (büyük sistemde yük olmasın: ≤ 65 bilinmeyen).
    K = [[float(v) for v in row] for row in result['K'].evalf(17).tolist()] if result['dofs'] <= 65 else None
    return {'text': text, 'tex': build_tex(text), 'nodal': [float(sp.sympify(v).evalf()) for v in result['a']],
            'exact': str(solution) if solution is not None else None,
            'exact_js': exact_for_browser(solution, result['x0'], result['x1']), 'K': K, **verdict_payload(result['exact_verdict']),
            'fixed': {str(i): float(v) for i, v in result['fixed'].items()},
            'first_order_residual': float(result['first_order']['residual']) if result['first_order'] else None,
            **classification_payload(result['classification'], result['degeneracy'])}


def exact_only(data):
    """Analitik çözüm + sürekli sınıflandırma (kâğıt yok): grafikler ve teşhis kutusu için hızlı istek. Ağ n = 1, p = 1
    ile en ucuz haliyle kurulur. Sürekli düzeyde kanıtlanmış ret (K1/K2) 'continuous_error' olarak döner; n = 1 ağının
    ayrık tekilliği sürekli problem hakkında bilgi vermediği için sınıflandırmayı silmez."""
    try:
        result = solve_problem({**data, 'n': 1, 'degree': 1}, symbolic=False, classify=True)
    except Diagnosis as error:
        if error.level == 'sürekli problem':
            return {'exact_js': None, 'continuous_error': {'category': error.category, 'message': str(error)}}
        return {'exact_js': None, **classification_payload(getattr(error, 'classification', None), getattr(error, 'degeneracy', None))}
    except UserError:
        return {'exact_js': None}
    verdict = exact_solution(result)
    return {'exact_js': exact_for_browser(verdict['u'], result['x0'], result['x1']), **verdict_payload(verdict),
            **classification_payload(result['classification'], result['degeneracy'])}


if __name__ == '__main__':
    payload = json.load(sys.stdin)
    try:
        if payload.get('mode') == 'exact':
            output, code = exact_only(payload), 0
        else:
            output, code = make_sheet(payload, payload.get('symbolic', True)), 0
    except UserError as error:
        output, code = {'error': str(error)}, 1
    sys.stdout.write(json.dumps(output, ensure_ascii=False))
    sys.stdout.flush()
    os._exit(code)      # süre sınırını aşmış bir dsolve iş parçacığı varsa beklemeden çık
