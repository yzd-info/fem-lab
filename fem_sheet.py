#!/usr/bin/env python3
"""fem_sheet.py — Galerkin sonlu eleman hattı ve çözüm kâğıdı (SymPy ile tam/sembolik hesap).

Problem:   (α(x) u')' + β(x) u' + γ(x) u = f(x),   x ∈ (x0, x1)
Çıktı:     stdin'den JSON alır, stdout'a JSON yazar: { text, tex, nodal, exact, exact_js, K, fixed, trace, verification, … }
           (mode = 'exact' ise yalnızca analitik çözüm ve sınıflandırma: grafiklerde bağımsız referans)
           text = "## başlık / $$formül$$" blokları (sayfada MathJax ile çizilir), tex = derlenebilir LaTeX belgesi.

Beş katman (veri yalnız aşağı akar; her katman bir öncekinin YAPILANDIRILMIŞ çıktısını alır):
  Katman 1  Problem tanımı       define_problem    kullanıcının metni olduğu gibi (tek kaynak; değiştirilmez)
  Katman 2  Problem doğrulaması  validate_problem  metin → iç gösterim + eşdeğerlik izi (iki bağımsız ayrıştırıcı)
  Katman 3  FEM motoru           fem_engine        şekil fonksiyonu, integraller, montaj, sınır koşulu, çözüm, yeniden kurulum
  Katman 4  Bağımsız doğrulama   fem_verify.verify ayrı dosya; motorun fonksiyonlarını ÇAĞIRMAZ, çıktısını bağımsız yollarla sınar
  Katman 5  Çözüm kâğıdı         build_sheet       yalnız 2–4'ün verisini biçimler; hiçbir sayı hesaplamaz
  run_pipeline bu sırayı uygular; bir katman başarısızsa sonrakiler çalışmaz ve neden Türkçe bildirilir.

Neden SymPy: integrallerin TEK TEK ve TAM alınması bir bilgisayar cebiri sistemi ister. Tarayıcıdaki JavaScript
çekirdeği yalnızca canlı grafikler için sayısal çözer; kâğıttaki her sayı bu hattan gelir.

Dosya haritası (yukarıdan aşağıya, her bölüm yalnızca üsttekileri kullanır):
  1. Girdi okuma      — metin → SymPy ifadesi (güvenli, beyaz listeli) ve bağımsız ikinci ayrıştırıcı
  2. Biçimleyiciler   — sayı / matris / denklem → LaTeX
  3. Katman 1–2       — problem tanımı, doğrulama, eşdeğerlik izi, katsayı dejenereliği
  4. Katman 3         — FEM motoru (yapılandırılmış veri üretir; LaTeX üretmez)
  5. Katman 5         — kâğıt (yalnız biçimleme)
  6. Hat, LaTeX belgesi ve komut satırı girişi (Katman 4 burada fem_verify'den çağrılır)
"""
import json
import os
import re
import sys
import threading

import mpmath
import sympy as sp
from sympy.parsing.sympy_parser import parse_expr, standard_transformations, implicit_multiplication_application, convert_xor

x, s = sp.symbols('x s', real=True)

# ---------- 1. GİRDİ OKUMA ----------
# Yalnızca bu adlara izin verilir; metin eval'den önce simge simge denetlendiği için başka hiçbir şey çalıştırılamaz.
NAMES = {
    'x': x, 'pi': sp.pi, 'e': sp.E, 'sin': sp.sin, 'cos': sp.cos, 'tan': sp.tan, 'asin': sp.asin, 'acos': sp.acos,
    'atan': sp.atan, 'sinh': sp.sinh, 'cosh': sp.cosh, 'tanh': sp.tanh, 'exp': sp.exp, 'ln': sp.log, 'log': sp.log,
    'sqrt': sp.sqrt, 'abs': sp.Abs,
}
DEGREE_NAMES = {1: 'linear', 2: 'quadratic', 3: 'cubic', 4: 'quartic'}
MATRIX_CAP = 9        # bundan büyük matrisler ilk 4 / son 2 satır-sütunla gösterilir
DETAIL_LIMIT = 8      # eleman sayısı bunu aşarsa integraller yalnızca ilk 2 ve son elemanda tek tek yazılır
EXACT_SOLVE_LIMIT = 8 # serbest bilinmeyen sayısı bunu aşarsa (ya da K rasyonel değilse) sistem ondalıklı çözülür (tam çözüm sayfalarca kesir üretir)


class UserError(ValueError):
    """Kullanıcının düzeltebileceği girdi hatası; mesajı aynen arayüzde gösterilir."""


# Teşhis kategorileri (tarayıcıdaki CATEGORY ile aynı adlar; iki çekirdek bağımsız karar verir). Tek bir genel "tekil"
# mesajı yerine NEDEN ve DÜZEY (sürekli problem / ayrık sistem) söylenir; ayrık tekillik, diferansiyel problemin
# çözümsüz olduğu anlamına gelmez ve öyle sunulmaz.
CATEGORY = {'K1': 'No solution (inconsistent system)', 'K2': 'Non-unique solution (null space ≠ {0})',
            'K3': 'Nearly singular / ill-conditioned', 'K4': 'Degenerate coefficient', 'K5': 'Numerical solver failure'}


# Teşhis düzeyi: makine değeri ('continuous' | 'discrete') → mesajda görünen. Neden ayrı: dil değişse de sunucu,
# testler ve korpus çalıştırıcısı düzeyi metinden ayrıştırmak zorunda kalmasın.
LEVEL_LABEL = {'continuous': 'continuous problem', 'discrete': 'discrete system'}
# Kanıt düzeyi / yöntem: makine değeri → ekranda görünen. Neden ayrı: bu değerler mantıkta karşılaştırılır (ör. yalnız
# sembolik kanıtlı K1/K2 ret); dil değişse de karşılaştırma bozulmasın. Tabloda olmayan (ölçülen değer gibi) metin olduğu gibi yazılır.
EVIDENCE_LABEL = {'symbolic': 'symbolic', 'numeric': 'numerical', 'numeric_50': 'numerical (50 digits)',
                  'numeric_9pt': 'numerical (9 points, 30 digits)', 'not_evaluable': 'could not be evaluated'}


class Diagnosis(UserError):
    """Kategorili teşhis; UserError olduğu için mevcut hata yolu (Türkçe mesaj, sunucu 400) aynen çalışır."""
    def __init__(self, category, level, detail):
        super().__init__(f'{CATEGORY[category]} — {LEVEL_LABEL[level]}: {detail}')
        self.category, self.level = category, level


# Sezgisel sayısal toleranslar (evrensel sınır DEĞİL; göreli tanımlı, tarayıcıdakilerle aynı değerler, testlerle sınandı).
RANK_TOL, COMPAT_TOL, COND_WARN, COND_SINGULAR = 1e-12, 1e-10, 1e10, 1e14
SVD_ANALYSIS_LIMIT = 65     # bundan büyük sayısal tekil sistemde mpmath SVD çok yavaş: analiz yapılmadığı açıkça söylenir
IVP_COARSE_DETAIL = ("the continuous problem has a unique solution (linear initial value problem), but the discrete system is singular or very badly "
                     "conditioned: the mesh is too coarse for this solution (for an oscillatory solution ωh must be small; ωh < 2√3 for p = 1). Increase n or p.")


def parse(text, label):
    """'2x^2 + sin(pi*x)' → SymPy. Tarayıcıdaki ayrıştırıcıyla aynı kurallar: ^ üs, örtük çarpım (2x), ondalıklar tam kesir."""
    text = (str(text) if text is not None else '').strip().lower() or '0'
    if re.search(r'\d\.?e[+-]?\d', text):        # "1e-8" = 1·e − 8 okunurdu (e Euler sayısı): sessizce yanlış katsayı
        raise UserError(f'{label}: scientific notation (e.g. 1e-8) is not supported; "e" is Euler\'s number. Write 10^-8 instead.')
    out, previous_is_operand = [], False
    for token in re.findall(r"\d+\.?\d*|\.\d+|[a-z]+|\S", text):
        is_number = token[0].isdigit() or (token[0] == '.' and len(token) > 1)
        is_name = token.isalpha()
        if is_name and token not in NAMES:
            raise UserError(f'{label}: unknown symbol "{token}"')
        if not (is_number or is_name or token in list('+-*/^()')):
            raise UserError(f'{label}: unexpected symbol "{token}"')
        if previous_is_operand and (is_number or is_name or token == '('):
            out.append('*')                                   # örtük çarpım: 2x, 3sin(x), (x+1)(x-1)
        out.append(f"Rational('{token}')" if is_number else '**' if token == '^' else token)
        previous_is_operand = is_number or token == ')' or token in ('x', 'pi', 'e')
    try:
        return sp.sympify(eval(''.join(out), {'__builtins__': {}}, {**NAMES, 'Rational': sp.Rational}))
    except Exception:
        raise UserError(f'{label}: the expression could not be read')


def parse_number(text, label):
    value = parse(text, label)
    if value.free_symbols or not value.is_finite:
        raise UserError(f'{label}: must be a constant number')
    return value


SECOND_GRAMMAR = standard_transformations + (implicit_multiplication_application, convert_xor)


# Katman 2'nin ikinci ayrıştırıcısı. Neden: "girdi = iç gösterim" iddiası, girdiyi üreten ayrıştırıcının kendisiyle
# denetlenemez. SymPy'nin kendi dilbilgisi (örtük çarpım, öncelik, ondalık) bizim simge simge çeviricimizden bağımsızdır;
# ad tablosu (NAMES) ortaktır, çünkü 'ln' ya da 'e'nin anlamı bir tanımdır, ayrıştırma değil. Yalnız beyaz liste
# denetiminden GEÇMİŞ metne uygulanır (parse önce çalışır), bu yüzden güvenlik denetimini gevşetmez.
def independent_parse(text):
    text = (str(text) if text is not None else '').strip().lower() or '0'
    return sp.sympify(parse_expr(text, local_dict=dict(NAMES), transformations=SECOND_GRAMMAR))


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


# ---------- 3. KATMAN 1–2: PROBLEM TANIMI VE DOĞRULAMASI ----------
COEFFICIENTS = (('alpha', 'α(x)'), ('beta', 'β(x)'), ('gamma', 'γ(x)'), ('f', 'f(x)'))
CONDITION_KINDS = ('none', 'dirichlet', 'neumann', 'flux', 'both')
EQUIVALENCE_TOL = 1e-12     # iki ayrıştırıcının sayısal karşılaştırması (göreli; ondalık girdide Float/kesir farkı ~1e-17)


# Katman 1. Kullanıcı problemin TEK kaynağıdır. Neden metin olarak saklanıyor: sonraki her katmanın "kullanıcının
# yazdığı" ile karşılaştırılabileceği değişmez bir kayıt gerekir. Hiçbir alan yorumlanmaz, kırpılmaz, düzeltilmez.
def define_problem(data):
    text = lambda v: '' if v is None else str(v)
    raw = {key: text(data.get(key)) for key in ('alpha', 'beta', 'gamma', 'f', 'x0', 'x1')}
    for side in ('left', 'right'):
        c = data.get(side) or {}
        raw[side] = {'type': text(c.get('type') or 'none'), 'value': text(c.get('value')), 'slope': text(c.get('slope'))}
    raw['n'], raw['degree'] = data.get('n', 4), data.get('degree', 1)
    return raw


# İki ifadenin aynı olduğunun kanıtı: önce sembolik fark (sıfırsa kesin), değilse aralığın 9 iç noktasında 30 basamakta
# göreli fark. Neden sayısal yol da var: ondalık girdide (0.05) ikinci ayrıştırıcı Float verir, sembolik fark tam sıfır olmaz.
def same_expression(first, second, x0, x1):
    try:
        if sp.simplify(first - second) == 0:
            return True, 'symbolic'
    except Exception:
        pass
    with mpmath.workdps(30):
        for k in range(9):
            pt = x0 + (x1 - x0) * sp.Rational(2 * k + 1, 18)
            try:
                a, b = complex(sp.N(first.subs(x, pt), 30)), complex(sp.N(second.subs(x, pt), 30))
            except Exception:
                return False, 'not_evaluable'
            if not abs(a - b) <= EQUIVALENCE_TOL * max(1.0, abs(a), abs(b)):
                return False, f'x = {dec(pt)}: {dec(a.real)} ≠ {dec(b.real)}'
    return True, 'numeric_9pt'


# Önde gelen katsayının (2. mertebede α, 1. mertebede β) sıfırları. Neden: α'nın sıfırı ya da işaret değişimi standart
# H¹ kuramını (düzgün eliptiklik |α| ≥ α₀ > 0) bozar; ama TEK BAŞINA problemin çözümsüz olduğunu göstermez. Bu yüzden sonuç
# yalnız sınıflandırma ve uyarıdır, ret değildir. Sıfırlar mümkünse sembolik (solveset), değilse sayısal bulunur ve
# hangi yolun kullanıldığı kaydedilir. Sayısal yol çok dar bir değme sıfırını kaçırabilir (bilinen sınır).
# Katman 2'de durur: problemin bir ÖZELLİĞİDİR (sınıfını belirler), çözümün varlığı hakkında hüküm vermez.
def coefficient_degeneracy(lead, x0, x1):
    span = x1 - x0
    if not lead.free_symbols:
        return {'kind': 'uniform', 'points': [], 'endpoints': [], 'interior': [], 'sign_change': [], 'method': 'symbolic'}
    points, method = None, 'numeric'
    try:
        if lead.is_rational_function(x):
            found = sp.solveset(sp.numer(sp.together(lead)), x, sp.Interval(x0, x1))
            if isinstance(found, sp.FiniteSet):
                points, method = sorted(found, key=float), 'symbolic'
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


# Katman 2. Amaç: çözücüye giden problemin kullanıcının yazdığıyla AYNI olduğunu göstermek. Çözümün varlığı/tekliği burada
# SORULMAZ (o Katman 4'ün işi). Her metin alanı iki bağımsız ayrıştırıcıdan geçer (parse: simge simge çevirici;
# independent_parse: SymPy dilbilgisi) ve sonuç iz tablosuna yazılır. Kural ihlali ya da eşdeğerliğin kurulamaması hesabı
# DURDURUR (UserError). Girdi asla sessizce düzeltilmez; büyük/küçük harf birleştirmesi tek dönüşümdür ve izde yazılır.
def validate_problem(raw):
    trace = []

    def row(item, given, internal, ok, how):
        trace.append({'item': item, 'raw': given, 'internal': internal, 'ok': bool(ok), 'how': how})
        if not ok:
            raise UserError(f'Problem validation: the input {item} and its internal representation are not equivalent ({how}). Computation stopped.')

    def second(text, label):
        try:
            return independent_parse(text)
        except Exception:
            raise UserError(f'Problem validation: {label} could not be read by the second parser; equivalence could not be established. Computation stopped.')

    shown = lambda text: text if text.strip() else ''
    x0, x1 = parse_number(raw['x0'], 'x₀'), parse_number(raw['x1'], 'x₁')
    for label, text, value in (('x₀', raw['x0'], x0), ('x₁', raw['x1'], x1)):
        ok, how = same_expression(value, second(text, label), 0, 1)
        row(label, shown(text), sp.latex(value), ok, how)
    if not x0 < x1:
        raise UserError('The interval requires x₀ < x₁.')

    coef = {}
    for key, label in COEFFICIENTS:
        coef[key] = parse(raw[key], label)
        ok, how = same_expression(coef[key], second(raw[key], label), x0, x1)
        row(label, shown(raw[key]), sp.latex(coef[key]), ok, how)
    alpha, beta, gamma, f = coef['alpha'], coef['beta'], coef['gamma'], coef['f']

    def whole(value):                         # 4.5 sessizce 4'e indirilmesin: tam sayı değilse durdur
        try:
            number = float(value)
        except (TypeError, ValueError):
            raise UserError('n and p must be integers.')
        if not number.is_integer():
            raise UserError('n and p must be integers.')
        return int(number)
    n, degree = whole(raw['n']), whole(raw['degree'])
    if not (1 <= n <= 64 and 1 <= degree <= 4):
        raise UserError('n must be in 1–64 and p in 1–4.')
    order = 2 if alpha != 0 else 1 if beta != 0 else 0
    if order == 0:
        raise UserError('α and β cannot both be zero: this is not a differential equation.')

    # Koşul türleri: none · dirichlet (u = g) · neumann (u' = g, TÜREV koşulu) · flux (α u' = q, AKI koşulu) ·
    # both (aynı uçta u ve u': başlangıç/son değer problemi; u' değeri slope alanında). Türev ve akı koşulu birbirine
    # dönüştürülmez: zayıf forma giren akı neumann'da Q = α(x_b)·g, flux'ta Q = q'dur ve iz tablosunda açıkça yazılır.
    conditions, slopes = {}, {}
    for side, point in (('left', x0), ('right', x1)):
        given = raw[side]
        kind = given['type']
        if kind not in CONDITION_KINDS:
            raise UserError('Unknown condition type.')
        value = parse_number(given['value'], 'condition value') if kind != 'none' else None
        conditions[side] = (kind, value)
        if kind == 'both':
            slopes[side] = parse_number(given['slope'], "u' condition value")
    kinds = [conditions['left'][0], conditions['right'][0]]
    if order == 1 and (any(k in kinds for k in ('neumann', 'flux', 'both')) or kinds.count('dirichlet') != 1):
        raise UserError('A first-order equation takes exactly one condition: give a value of u at one end and leave the other end empty.')
    # Başlangıç değer problemi: iki koşul aynı uçta. Öbür uca koşul verilirse problem aşırı belirlenir (3 koşul, 2. mertebe).
    ivp = 'both' in kinds
    if ivp and sorted(kinds) != ['both', 'none']:
        raise UserError("If both u and u' are given at one end (initial value problem), no condition may be given at the other end.")
    # u' (ya da akı) koşulu zayıf forma yalnızca akı Q = α u' üzerinden girer. α(uç) = 0 ise verilen değer sessizce
    # kaybolurdu (u' koşulunda) ya da 0 = q çelişkisi olurdu (akı koşulunda). Bu yüzden o uçta α ≠ 0 şartı aranır.
    for side, point in (('left', x0), ('right', x1)):
        if conditions[side][0] in ('neumann', 'both') and alpha.subs(x, point) == 0:
            raise UserError("α is zero at this end: a u' condition does not enter the equation (the flux α u' vanishes).")
        if conditions[side][0] == 'flux' and alpha.subs(x, point) == 0:
            raise UserError("α is zero at this end: a flux condition α u' = q is meaningful only for q = 0, and then it is the natural condition itself; "
                            "leave this end as 'not given' or check α.")

    for side, point, name in (('left', x0, 'condition at x₀'), ('right', x1, 'condition at x₁')):
        kind, value = conditions[side]
        given, P = raw[side], sp.latex(point)
        a_b = alpha.subs(x, point)
        if kind == 'none':
            meaning = ('\\text{no condition (first order)}' if order == 1 else '\\text{free end (initial value problem)}' if ivp
                       else '\\alpha(' + P + ') = 0\\ \\text{(degenerate end; the flux vanishes by itself)}' if a_b == 0
                       else f"\\alpha u'({P}) = 0\\ \\text{{(natural condition)}}")
            trace.append({'item': name, 'raw': '(not given)', 'internal': meaning, 'ok': True, 'how': 'type: not given'})
            continue
        ok, how = same_expression(value, second(given['value'], name), x0, x1)
        if kind == 'dirichlet':
            internal = f'u({P}) = {sp.latex(value)}\\ \\text{{(essential; flux unknown)}}'
        elif kind == 'neumann':
            internal = f"u'({P}) = {sp.latex(value)}\\ \\Rightarrow\\ Q = \\alpha({P})\\,u'({P}) = {sp.latex(sp.simplify(a_b * value))}"
        elif kind == 'flux':
            internal = f"\\alpha({P})\\,u'({P}) = {sp.latex(value)}\\ \\Rightarrow\\ Q = {sp.latex(value)}"
        else:
            ok2, how2 = same_expression(slopes[side], second(given['slope'], name), x0, x1)
            ok, how = ok and ok2, how if ok2 else how2
            internal = (f"u({P}) = {sp.latex(value)},\\ u'({P}) = {sp.latex(slopes[side])}\\ \\Rightarrow\\ "
                        f"Q = {sp.latex(sp.simplify(a_b * slopes[side]))}")
        label = {'dirichlet': 'u =', 'neumann': "u' =", 'flux': "α u' =", 'both': "u and u' ="}[kind]
        row(name, f"{label} {given['value']}" + (f", {given['slope']}" if kind == 'both' else ''), internal, ok, how)

    present = ', '.join(f'{name} {"≢" if coef[key] != 0 else "≡"} 0' for key, name in (('alpha', 'α'), ('beta', 'β'), ('gamma', 'γ'), ('f', 'f')))
    trace.append({'item': 'operator', 'raw': f"(({shown(raw['alpha']) or '0'}) u')' + ({shown(raw['beta']) or '0'}) u' + "
                  f"({shown(raw['gamma']) or '0'}) u = {shown(raw['f']) or '0'}",
                  'internal': tex_equation(alpha, beta, gamma, f), 'ok': True, 'how': f'{present}; order {order}'})
    trace.append({'item': 'mesh', 'raw': f"n = {raw['n']}, p = {raw['degree']}", 'internal': f'n = {n},\\ p = {degree}',
                  'ok': True, 'how': 'integer'})

    degeneracy = coefficient_degeneracy(alpha if order == 2 else beta, x0, x1)
    problem = dict(alpha=alpha, beta=beta, gamma=gamma, f=f, x0=x0, x1=x1, n=n, degree=degree, order=order,
                   conditions=conditions, slopes=slopes, ivp=ivp, degeneracy=degeneracy)
    return {'definition': raw, 'problem': problem, 'trace': trace}


# ---------- 4. KATMAN 3: FEM MOTORU ----------
# Motor yalnız YAPILANDIRILMIŞ veri üretir (SymPy nesneleri, sözlükler); LaTeX metni üretmez. Kâğıt (Katman 5) bu veriyi
# biçimler, doğrulama (Katman 4, fem_verify.py) aynı veriyi motorun kullanmadığı yollarla sınar.
GLOBAL_SHOW = 2       # eleman başına ayrıntı (global N, u_h, R) en fazla bu kadar elemanda; diğerleri aynı yapıda
ROW_SHOW = 9          # test fonksiyonu başına denklem bu kadar satıra kadar tam yazılır


def lagrange_basis(degree, h):
    """Yerel koordinat s ∈ [0, h] üzerinde eşit aralıklı düğümlü Lagrange şekil fonksiyonları ve türevleri."""
    nodes = [sp.Rational(k, degree) * h for k in range(degree + 1)]
    N = [sp.expand(sp.prod([(s - nodes[m]) / (nodes[k] - nodes[m]) for m in range(degree + 1) if m != k]))
         for k in range(degree + 1)]
    return N, [sp.diff(Nk, s) for Nk in N]


def hat(k):
    return sp.Symbol(f'a^{{(e)}}_{{{k}}}')      # e elemanının yerel katsayısı a_k^{(e)} (global a_j ile karıştırılmasın)


def node_value(j):
    return sp.Symbol(f'a_{{{j}}}')


# Şekil fonksiyonlarının polinom interpolasyonundan kurulması (Vandermonde yolu): u^e(s) = Σ c_m s^m, u^e(s_k) = â_k.
# Neden motorda: N_k'nın NEREDEN geldiği hesabın parçasıdır. Ürün formülüyle (lagrange_basis) kurulan N ile karşılaştırılır;
# Kronecker-delta, birim bölüşüm ve polinom tamlığı (Σ_k N_k(s) s_k^m = s^m, m ≤ p) burada sınanır. Uyuşmazlık hesabı durdurur.
def shape_function_construction(degree, h, N_engine):
    cs = [sp.Symbol(f'c_{{{m}}}') for m in range(degree + 1)]
    hats = [hat(k + 1) for k in range(degree + 1)]
    nodes = [sp.Rational(k, degree) * h for k in range(degree + 1)]
    poly = sum(c * s ** m for m, c in enumerate(cs))
    conds = [sp.expand(poly.subs(s, v)) for v in nodes]
    V = sp.Matrix([[v ** m for m in range(degree + 1)] for v in nodes])
    sol = sp.solve([sp.Eq(cond, hv) for cond, hv in zip(conds, hats)], cs, dict=True)[0]
    shapes = [sp.expand(sp.expand(poly.subs(sol)).coeff(hv)) for hv in hats]
    for k, (Nk, Ne) in enumerate(zip(shapes, N_engine)):
        if sp.simplify(Nk - Ne) != 0:
            raise RuntimeError(f'The derived N_{k + 1} differs from the one used in the computation')
    delta = sp.Matrix(degree + 1, degree + 1, lambda i, j: shapes[i].subs(s, nodes[j]))
    partition = sp.simplify(sum(shapes))
    completeness = [sp.simplify(sum(Nk * nodes[k] ** m for k, Nk in enumerate(shapes)) - s ** m) == 0 for m in range(degree + 1)]
    if delta != sp.eye(degree + 1) or partition != 1 or not all(completeness):
        raise RuntimeError('The shape functions violate δ_ij / partition of unity / completeness')
    return dict(degree=degree, h=h, cs=cs, hats=hats, nodes=nodes, poly=poly, conds=conds, V=V, sol=sol,
                W=V[1:, 1:], rhs=[hv - hats[0] for hv in hats[1:]], shapes=shapes,
                derivatives=[sp.diff(Nk, s) for Nk in shapes], delta=delta, partition=partition, completeness=completeness)


# Tek bir eleman integrali. Dönüş YAPILANDIRILMIŞ kayıttır: çarpanlar, açılmış integrand, ilkel fonksiyon, değer, yöntem
# ('zero' / 'polynomial' / 'symbolic' / 'numeric'). Kapalı biçim bulunamazsa mpmath.quad (uyarlamalı tanh–sinh) ile 15 basamak.
# Önbellek: aynı integrand (SymPy kanonik biçimi; çarpım değişmeli olduğundan β = 0 iken N_i'αN_j' ile N_j'αN_i' AYNI
# ifadedir) aynı h ve yolla aynı sonucu verir; yeniden hesaplanmaz, kayıt kendi etiketiyle tek tek yazılır. Neden: aşkın
# katsayıda (ör. α = √(1+x)) her integral saniyeler sürer; simetrik K'de yarısı tekrardı (korpus C32 sunucu sınırını aşıyordu).
_INTEGRAL_CACHE = {}


def integral(label, terms, h, symbolic):
    integrand = sum(sign * sp.prod(factors) for sign, factors in terms)
    key = (integrand, h, bool(symbolic))
    if key in _INTEGRAL_CACHE:
        return {**_INTEGRAL_CACHE[key], 'label': label, 'terms': terms}
    result = _integral(label, terms, h, symbolic, integrand)
    _INTEGRAL_CACHE[key] = result
    return dict(result)                         # kopya: çağıranın kaydı değiştirmesi önbelleği bozmasın


# integral()'in hesabı (önbelleksiz). Ayrı fonksiyon: önbellek anahtarı ile hesap birbirine karışmasın.
def _integral(label, terms, h, symbolic, integrand):
    record = {'label': label, 'terms': terms, 'h': h, 'integrand': None, 'primitive': None}
    polynomial = integrand.is_polynomial(s)
    if polynomial:
        integrand = sp.expand(integrand)
    if integrand == 0:
        return {**record, 'method': 'zero', 'value': sp.Integer(0)}
    if symbolic:
        try:
            if polynomial:
                primitive = sp.integrate(integrand, s)
                value = sp.simplify(primitive.subs(s, h) - primitive.subs(s, 0))
                return {**record, 'method': 'polynomial', 'integrand': integrand, 'primitive': primitive, 'value': value}
            value = sp.simplify(sp.integrate(integrand, (s, 0, h)))
            if not value.has(sp.Integral) and value.is_finite and value.is_real:
                # İlkel fonksiyon yalnız gösterim içindir: F(h) − F(0) değere sembolik olarak eşit değilse kâğıda yazılmaz.
                # Açılmamış integrand integre edilir (aynı fonksiyon; açılmış biçimin belirsiz integrali √ ve trigonometrik
                # çarpımlarda birkaç kat yavaştı); kâğıtta açılmış integrand gösterilir, eşitlik aşağıdaki denetimle güvencededir.
                expanded = sp.expand(integrand)
                try:
                    primitive = sp.integrate(integrand, s)
                    ok = not primitive.has(sp.Integral) and sp.simplify(primitive.subs(s, h) - primitive.subs(s, 0) - value) == 0
                except Exception:
                    ok = False
                return {**record, 'method': 'symbolic', 'value': value,
                        'integrand': expanded if ok else None, 'primitive': primitive if ok else None}
        except Exception:
            pass
    value = sp.Float(mpmath.quad(sp.lambdify(s, integrand, 'mpmath'), [0, float(h)]), 15)
    return {**record, 'method': 'numeric', 'value': value}


# Tam (rasyonel/sembolik) yolda tekil sistemin sınıflandırılması: rank ve uyumluluk KESİN (tolerans yok).
def classify_discrete_exact(K_ff, F_f, ivp, has_dirichlet):
    if ivp:
        return Diagnosis('K3', 'discrete', IVP_COARSE_DETAIL)
    rank, rank_aug, size = K_ff.rank(simplify=True), K_ff.row_join(F_f).rank(simplify=True), K_ff.shape[1]
    if rank_aug > rank:
        example = '' if has_dirichlet else " E.g. pure Neumann with ∫f ≠ flux difference."     # örnek yalnız u hiç verilmemişse ilgili
        return Diagnosis('K1', 'discrete', f'rank K = {rank} < rank [K|F] = {rank_aug} (exact arithmetic): the right-hand side is inconsistent. '
                                                f'This is a verdict on the DISCRETE system; for the continuous problem see the classification.{example}')
    null = K_ff.nullspace()[0]
    shape = 'constant vector (u + c is also a solution)' if len(set(sp.simplify(v / null[0]) for v in null)) == 1 else 'a non-constant vector'
    advice = '' if has_dirichlet else " Prescribing u at one end or changing γ may remove the null space."   # öneri yalnız u hiç verilmemişse anlamlı
    return Diagnosis('K2', 'discrete', f'rank K = rank [K|F] = {rank} < {size} (exact arithmetic): null-space dimension {size - rank}, '
                                            f'example: {shape}.{advice}')


# Sayısal yolda tekil/neredeyse tekil sistem: mpmath SVD ile rank (σ ≤ RANK_TOL·σ_max) ve normwise geriye hata
# η = ‖b − Ax̂‖/(‖A‖‖x̂‖ + ‖b‖) ile uyumluluk. Büyük sistemde analiz yapılmadığı açıkça söylenir (iddia yok).
def classify_discrete_numeric(A, b, ivp, condition):
    if ivp:
        return Diagnosis('K3', 'discrete', IVP_COARSE_DETAIL)
    if A.rows > SVD_ANALYSIS_LIMIT:
        return Diagnosis('K3', 'discrete', f'numerically singular (κ₁ ≈ {condition:.1e}); rank/null-space analysis not performed for {A.rows} unknowns.')
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
            return Diagnosis('K2', 'discrete', f'rank {rank} < {A.cols}; the right-hand side is consistent (η = {eta:.1e}), null-space dimension {A.cols - rank}.')
        return Diagnosis('K1', 'discrete', f'rank {rank} < {A.cols}; the right-hand side is not in the range (η = {eta:.1e} > {COMPAT_TOL}).')
    return Diagnosis('K3', 'discrete', f'numerically singular: full rank but κ₁ ≈ {condition:.1e}; the result is unreliable.')


# Bant matris LU'su. Neden var: indirgenmiş FEM matrisi bantlıdır (eleman başına p + 1 düğüm, yarı bant ≤ p + 1), ama
# mpmath'in yoğun LU_decomp'u ve cond'u (ters matrisle) O(n³) saf Python'dur: n = 64, p = 4'te (255 bilinmeyen) ~40 s sürüyor,
# sunucunun 60 s sınırını aşıp HTTP 504 veriyordu (korpus C48/C50). Algoritma mpmath.LU_decomp'un AYNISIDIR — ölçekli kısmi
# pivotlama (|A_kj| / Σ_l |A_kl|, eşitlikte ilk satır), aynı tekillik eşiği |pivot| ≤ ‖A‖₁·ε → ZeroDivisionError — yalnız
# bant dışındaki kesin sıfırlar üzerinde döngü kurulmaz. Pivotlamada üst bant en çok alt bant kadar genişler (LAPACK gbtrf),
# bu yüzden sınırlar kl ve ku + kl'dir. Maliyet O(n·kl·(ku + kl)). Satırlar mpf listeleri: mpmath.matrix girdi erişimi yavaştır.
def band_lu(rows):
    n = len(rows)
    A = [list(row) for row in rows]
    kl = max([i - j for i in range(n) for j in range(i) if A[i][j]] or [0])
    ku = max([j - i for i in range(n) for j in range(i + 1, n) if A[i][j]] or [0])
    width = ku + kl
    tol = mpmath.mnorm(mpmath.matrix(A), 1) * mpmath.eps if n else 0
    perm, first = [], list(range(n))           # first[i]: i. satırda saklı çarpanların başladığı sütun (ileri çözüm sınırı)
    for j in range(n - 1):
        last = min(n - 1, j + width)
        best, pivot = 0, None
        for k in range(j, min(n - 1, j + kl) + 1):
            s = mpmath.fsum(abs(A[k][l]) for l in range(j, last + 1))
            if s <= tol:
                raise ZeroDivisionError('matrix is numerically singular')
            current = abs(A[k][j]) / s
            if current > best:
                best, pivot = current, k
        if pivot is None:
            raise ZeroDivisionError('matrix is numerically singular')
        A[j], A[pivot] = A[pivot], A[j]
        first[j], first[pivot] = first[pivot], first[j]
        perm.append(pivot)
        if abs(A[j][j]) <= tol:
            raise ZeroDivisionError('matrix is numerically singular')
        for i in range(j + 1, min(n - 1, j + kl) + 1):
            m = A[i][j] = A[i][j] / A[j][j]
            first[i] = min(first[i], j)
            if m:
                Aj, Ai = A[j], A[i]
                for c in range(j + 1, last + 1):
                    Ai[c] -= m * Aj[c]
    if n and abs(A[n - 1][n - 1]) <= tol:
        raise ZeroDivisionError('matrix is numerically singular')
    return {'A': A, 'perm': perm, 'first': first, 'width': width, 'n': n}


# Bant LU ile A x = b (mpmath.L_solve + U_solve sırası: önce satır değişimleri, sonra birim alt üçgen, sonra üst üçgen).
# start: b'nin ilk sıfır olmayan (permütasyondan sonraki) konumu; birim vektör çözümlerinde ileri çözüm oradan başlar.
def band_solve(lu, b, start=0):
    A, n, width = lu['A'], lu['n'], lu['width']
    y = list(b)
    for j, k in enumerate(lu['perm']):
        y[j], y[k] = y[k], y[j]
    for i in range(start + 1, n):
        row = A[i]
        y[i] -= mpmath.fsum(row[j] * y[j] for j in range(max(lu['first'][i], start), i))
    for i in range(n - 1, -1, -1):
        row = A[i]
        y[i] = (y[i] - mpmath.fsum(row[j] * y[j] for j in range(i + 1, min(n, i + width + 1)))) / row[i]
    return y


# mpmath.lu_solve ile aynı sözleşme (+10 bit çalışma duyarlığı, tekilde ZeroDivisionError), bant yapısıyla.
def band_lu_solve(rows, b):
    with mpmath.extraprec(10):
        return band_solve(band_lu(rows), b)


# κ₁ = ‖A‖₁ ‖A⁻¹‖₁ (mpmath.cond'un varsayılanı, tahmin değil kesin değer). ‖A⁻¹‖₁ = max_j ‖A⁻¹ e_j‖₁: n bant çözümü,
# ters matrisin tamamı yazılmaz. e_j permütasyondan sonra tek bir konumda 1'dir; ileri çözüm o konumdan başlar.
def band_condition(rows):
    n = len(rows)
    if not n:
        return mpmath.mpf(1)
    with mpmath.extraprec(20):                      # mpmath.inverse'ün çalışma duyarlığı
        lu = band_lu(rows)
        where = list(range(n))
        for j, k in enumerate(lu['perm']):
            where[j], where[k] = where[k], where[j]
        position = {v: i for i, v in enumerate(where)}
        inv_norm = 0
        for j in range(n):
            e = [mpmath.mpf(0)] * n
            e[j] = mpmath.mpf(1)
            col = band_solve(lu, e, start=position[j])
            inv_norm = max(inv_norm, mpmath.fsum(abs(v) for v in col))
    return mpmath.mnorm(mpmath.matrix(rows), 1) * inv_norm


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
        z = band_lu_solve(A2T.tolist(), [mpmath.mpf(float(v)) for v in A.row(d)])   # bant: 1. mertebe n = 64, p = 4'te de saniyeler
        y = [mpmath.mpf(1) if i == d else -z[rows2.index(i)] for i in range(K.rows)]
        bn = [mpmath.mpf(float(v)) for v in b]
        yb = sum(yi * bi for yi, bi in zip(y, bn))
        scaled = float(abs(yb) / (mpmath.sqrt(sum(v * v for v in y)) * mpmath.sqrt(sum(v * v for v in bn))))
        compatible = scaled <= COMPAT_TOL
    return {'row': d + 1, 'residual': residual, 'yb': yb, 'scaled': scaled, 'compatible': compatible, 'exact': exact}


# İndirgenmiş sistemin çözümü (önceki gövde aynen; yalnız ayrı fonksiyona alındı ki ayrık teşhis hatası sınıflandırma
# bilgisiyle birlikte taşınabilsin). Tam yolda kesirli LU, sayısal yolda bant LU + κ₁ (mpmath.lu_solve / cond ile aynı
# aritmetik; tests.py ikisini karşılaştırır).
def solve_reduced(K_ff, F_f, exact, ivp_unique, has_dirichlet):
    if exact:
        if sp.simplify(K_ff.det()) == 0:
            raise classify_discrete_exact(K_ff, F_f, ivp_unique, has_dirichlet)
        return [sp.simplify(v) for v in K_ff.LUsolve(F_f)], None
    A = mpmath.matrix(K_ff.evalf(30).tolist())
    bvec = mpmath.matrix(F_f.evalf(30).tolist())
    rows = A.tolist()
    try:                                    # tam tekil matriste ZeroDivisionError: o da sınıflandırılır
        condition = float(band_condition(rows)) if K_ff.shape[0] else 1.0
    except ZeroDivisionError:
        condition = float('inf')
    if not condition <= COND_SINGULAR:
        raise classify_discrete_numeric(A, bvec, ivp_unique, condition)
    return [sp.Float(v, 15) for v in band_lu_solve(rows, [bvec[i] for i in range(bvec.rows)])], condition


# Elemanın global koordinattaki yaklaşımı: x = x_e + s, â_k = a_{(e−1)p+k}.
def element_approximation(r, e):
    degree, h, x0 = r['degree'], r['h'], r['x0']
    xe = x0 + e * h
    Ns = [sp.expand(Nk.subs(s, x - xe)) for Nk in r['N']]
    nodes = [node_value(e * degree + k + 1) for k in range(degree + 1)]
    return xe, Ns, nodes, sum(a * Nk for a, Nk in zip(nodes, Ns))


# Eleman İÇİNDEKİ klasik kalıntı: R_e = L[u_h^{(e)}] − f, klasik türevle; 2. mertebede (αu_h')' terimi de dahil. Global
# kalıntı bu değildir: C⁰ elemanda u_h' düğümlerde sıçrar ve (αu_h')' orada δ terimleri taşır (kâğıt ikisini ayırır).
def element_residual(r, uh):
    if r['order'] == 1:
        R = r['beta'] * sp.diff(uh, x) + r['gamma'] * uh - r['f']
    else:
        R = sp.diff(r['alpha'] * sp.diff(uh, x), x) + r['beta'] * sp.diff(uh, x) + r['gamma'] * uh - r['f']
    syms = sorted(uh.free_symbols - {x}, key=str)
    return sp.collect(sp.expand(R), syms)


# Katman 3'ün girişi: DOĞRULANMIŞ problem (validate_problem çıktısı). Motor girdiyi yeniden yorumlamaz; katsayılar,
# koşullar ve akı kuralı (neumann: Q = α(x_b)·g, flux: Q = q) aynen kullanılır.
def fem_engine(valid, symbolic=True):
    P = valid['problem']
    alpha, beta, gamma, f = P['alpha'], P['beta'], P['gamma'], P['f']
    x0, x1, n, degree, order = P['x0'], P['x1'], P['n'], P['degree'], P['order']
    conditions, slopes, ivp, degeneracy = P['conditions'], P['slopes'], P['ivp'], P['degeneracy']
    ivp_unique = ivp and degeneracy['kind'] == 'uniform'     # varlık-teklik teoremi α ≠ 0 iken geçerli

    h = (x1 - x0) / n
    dofs = n * degree + 1
    N, dN = lagrange_basis(degree, h)
    shape = shape_function_construction(degree, h, N)
    size = degree + 1
    uniform_K = not (alpha.free_symbols or beta.free_symbols or gamma.free_symbols)
    uniform_F = not f.free_symbols
    # Birinci mertebede (α = 0) kısmi integrasyon yoktur: ∫ N_i (β u' + γ u) dx = ∫ N_i f dx doğrudan yazılır.
    # İkinci mertebe formülün (−1) katıdır; çözüm aynıdır ama kâğıt ders kitabındaki işaretle ve doğru formülasyonla yazılır.
    sign = 1 if order == 2 else -1
    detailed = list(range(n)) if n <= DETAIL_LIMIT else [0, 1, n - 1]

    elements = []                                   # eleman başına: K^e, F_l^e ve (varsa) tek tek integral kayıtları
    for e in range(n):
        xe = x0 + e * h
        local = lambda g: g.subs(x, xe + s)         # katsayıyı elemanın yerel koordinatına taşı: x = x_e + s
        show = symbolic and e in detailed
        record = {'e': e, 'xe': xe, 'K_steps': [], 'F_steps': []}
        if uniform_K and e > 0:
            record['K'] = elements[0]['K']          # sabit katsayı + eşit eleman boyu → matris her elemanda aynı
        else:
            record['K'] = sp.zeros(size, size)
            for i in range(size):
                for j in range(size):
                    terms = ([(1, [dN[i], local(alpha), dN[j]]), (-1, [N[i], local(beta), dN[j]]), (-1, [N[i], local(gamma), N[j]])]
                             if order == 2 else [(1, [N[i], local(beta), dN[j]]), (1, [N[i], local(gamma), N[j]])])
                    step = integral(f'K^{{({e + 1})}}_{{{i + 1}{j + 1}}}', terms, h, show)
                    record['K'][i, j] = step['value']
                    record['K_steps'].append(step)
        if uniform_F and e > 0:
            record['F'] = elements[0]['F']
        else:
            record['F'] = sp.zeros(size, 1)
            for i in range(size):
                step = integral(f'F^{{({e + 1})}}_{{l,{i + 1}}}', [(-sign, [N[i], local(f)])], h, show)
                record['F'][i] = step['value']
                record['F_steps'].append(step)
        elements.append(record)

    K, F_l = sp.zeros(dofs, dofs), sp.zeros(dofs, 1)
    for record in elements:                         # montaj: yerel düğüm k → global e·p + k
        base = record['e'] * degree
        for i in range(size):
            F_l[base + i] += record['F'][i]
            for j in range(size):
                K[base + i, base + j] += record['K'][i, j]

    # Sınır terimi F_b = [N_i α u']: yalnızca uç düğümlerde sıfırdan farklıdır. Q0 = (αu')(x0), Q1 = (αu')(x1).
    # neumann (u' = g): Q = α(x_b)·g · flux (α u' = q): Q = q · both: Q = α(x_b)·u'(x_b) · none (BVP): Q = 0 (doğal koşul).
    # Q None = bu uçtaki akı bilinmiyor (u verilmiş uç ya da başlangıç değer probleminde koşulsuz kalan uç).
    Q = {}
    if order == 2:
        for side, point in (('left', x0), ('right', x1)):
            kind, value = conditions[side]
            Q[side] = (alpha.subs(x, point) * value if kind == 'neumann' else value if kind == 'flux'
                       else alpha.subs(x, point) * slopes[side] if kind == 'both'
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
    # Tam (kesirli) yol yalnız K_ff RASYONELKEN: o zaman det K tam rasyoneldir (tekillik kararı kesin ve anlık) ve çözüm
    # bileşenleri F'nin girdilerinin rasyonel katsayılı birleşimidir. K'de log/kök/π gibi aşkın girdiler varsa (ör. α = 1/(1+x)
    # → log) simplify(det K) ve çözümün sadeleştirilmesi dakikalarca sürüyordu (korpus C31: > 10 dk; sunucu 60 s sonra
    # sayısal yola düşüyordu) ve sonuç yüzlerce terimlik, okunmaz ifadeydi. O durumda integraller yine tam yazılır, yalnız
    # indirgenmiş sistem ondalıklı çözülür; neden solve_note'ta taşınır ve kâğıtta a'nın yanına yazılır.
    decimal_input = K_ff.has(sp.Float) or F_f.has(sp.Float)
    solve_note = ('decimal coefficient' if decimal_input else f'{len(free)} > {EXACT_SOLVE_LIMIT} free unknowns' if len(free) > EXACT_SOLVE_LIMIT
                  else 'K is not rational' if not all(v.is_Rational for v in K_ff) else None)
    exact = solve_note is None
    a_free, condition = solve_reduced(K_ff, F_f, exact, ivp_unique, bool(fixed))
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
    reaction_slopes = {side: sp.simplify(flux / alpha.subs(x, point)) for side, point, flux in
                       (('left', x0, reactions.get('left')), ('right', x1, reactions.get('right')))
                       if flux is not None and alpha.subs(x, point) != 0}           # u_h'(x_b) = Q / α(x_b)

    r = dict(alpha=alpha, beta=beta, gamma=gamma, f=f, x0=x0, x1=x1, n=n, degree=degree, order=order, h=h, dofs=dofs, first_order=first_order,
             conditions=conditions, slopes=slopes, ivp=ivp, N=N, dN=dN, shape=shape, elements=elements, detailed=detailed, uniform_K=uniform_K,
             uniform_F=uniform_F, K=K, F_l=F_l, F_b=F_b, Q=Q, fixed=fixed, free=free, rows=rows, K_ff=K_ff, F_f=F_f, a=a,
             exact=exact, solve_note=solve_note, reactions=reactions, reaction_slopes=reaction_slopes, symbolic=symbolic, condition=condition, degeneracy=degeneracy, classification=None)
    r['formulation'] = formulation_data(r)
    return r


# Galerkin formülasyonunun açık hali ve montaj / sınır koşulu dönüşümlerinin kaydı. Neden motorda: kâğıt bu ifadeleri
# YAZAR ama hesaplamaz. 1. mertebede ∫φ_i R dx doğrudan alınır ve K'nin satırıyla karşılaştırılır; ayrışırsa hesap durur
# (gösterilen türetme ile yapılan hesap ayrışamaz).
def formulation_data(r):
    degree, n, h, x0, dofs, order = r['degree'], r['n'], r['h'], r['x0'], r['dofs'], r['order']
    nodes = [x0 + sp.Rational(i, degree) * h for i in range(dofs)]
    a_syms = [node_value(j + 1) for j in range(dofs)]
    shown = r['detailed'][:GLOBAL_SHOW]
    with_rows = order == 1 and dofs <= ROW_SHOW
    elements = {}
    for e in sorted(set(shown) | (set(range(n)) if with_rows else set())):
        xe, Ns, coefficients, uh = element_approximation(r, e)
        elements[e] = {'xe': xe, 'N': Ns, 'coefficients': coefficients, 'uh': uh, 'residual': element_residual(r, uh)}

    galerkin_rows = []
    if with_rows:
        for i in r['rows']:
            parts, total = [], 0
            for e in range(n):
                k = i - e * degree
                if not 0 <= k <= degree:
                    continue
                E = elements[e]
                parts.append({'e': e, 'xe': E['xe'], 'weight': E['N'][k], 'residual': E['residual']})
                total += sp.integrate(sp.expand(E['N'][k] * E['residual']), (x, E['xe'], E['xe'] + h))
            total = sp.collect(sp.expand(total), a_syms)
            row = sp.expand(sum(r['K'][i, j] * a_syms[j] for j in range(dofs)) - r['F_l'][i])
            if total.has(sp.Integral) or (not r['K'].has(sp.Float) and sp.simplify(total - row) != 0):
                raise RuntimeError(f'Galerkin equation i={i + 1} does not match row i of K')
            galerkin_rows.append({'i': i, 'parts': parts, 'total': total})

    contributions, load_contributions, row_equations = [], [], []
    if dofs <= ROW_SHOW:
        shared = {}
        for record in r['elements']:
            base = record['e'] * degree
            for i in range(degree + 1):
                for j in range(degree + 1):
                    shared.setdefault((base + i, base + j), []).append((record['e'], i, j, record['K'][i, j]))
        contributions = [{'I': I, 'J': J, 'parts': c, 'value': r['K'][I, J]} for (I, J), c in sorted(shared.items()) if len(c) > 1]
        fshared = {}
        for record in r['elements']:
            for i in range(degree + 1):
                fshared.setdefault(record['e'] * degree + i, []).append((record['e'], i, record['F'][i]))
        load_contributions = [{'I': I, 'parts': c, 'value': r['F_l'][I]} for I, c in sorted(fshared.items()) if len(c) > 1]
        if order == 2:
            unknown = {0: -sp.Symbol('Q_0'), dofs - 1: sp.Symbol('Q_1')}
            for i in r['rows']:
                side = 'left' if i == 0 else 'right' if i == dofs - 1 else None
                Fb = r['F_b'][i] if side is None or r['Q'][side] is not None else unknown[i]
                row_equations.append({'i': i, 'lhs': sp.expand(sum(r['K'][i, j] * a_syms[j] for j in range(dofs))),
                                      'rhs': sp.expand(r['F_l'][i] + Fb)})

    reconstruction = []
    for e in r['detailed']:
        xe = x0 + e * h
        values = [r['a'][e * degree + k] for k in range(degree + 1)]
        piece = sp.expand(sum(v * Nk for v, Nk in zip(values, r['N'])).subs(s, x - xe))
        local = [sp.expand(Nk.subs(s, x - xe)) for Nk in r['N']] if e in shown else None
        reconstruction.append({'e': e, 'xe': xe, 'values': values, 'N': local, 'piece': piece})

    ivp_slope = None
    if r['ivp']:
        side = 'left' if r['conditions']['left'][0] == 'both' else 'right'
        e, s_at, point = (0, 0, x0) if side == 'left' else (n - 1, h, r['x1'])
        slope = sp.simplify(sum(sp.sympify(r['a'][e * degree + k]) * r['dN'][k].subs(s, s_at) for k in range(degree + 1)))
        ivp_slope = {'point': point, 'computed': slope, 'given': r['slopes'][side]}

    return dict(nodes=nodes, a_symbols=a_syms, shown=shown, elements=elements, galerkin_rows=galerkin_rows,
                contributions=contributions, load_contributions=load_contributions, row_equations=row_equations,
                reconstruction=reconstruction, ivp_slope=ivp_slope)


# ---------- 5. KATMAN 5: ÇÖZÜM KÂĞIDI ----------
# Kâğıt yalnız Katman 2 (iz tablosu, doğrulanmış problem), Katman 3 (motor kaydı) ve Katman 4 (doğrulama kaydı) verisini
# BİÇİMLER. Hiçbir sayı burada hesaplanmaz; sp.latex / tex yalnız gösterimdir. Bir değer yoksa "hesaplanmadı" yazılır.
STATUS_MARK = {'passed': '\\checkmark', 'failed': '\\times', 'warning': '!', 'not_run': '-'}
# Makine durum değeri (fem_verify) → kâğıtta görünen etiket. Neden ayrı: arayüz dili değişse de doğrulama verisi değişmesin.
STATUS_LABEL = {'passed': 'Passed', 'warning': 'Warning', 'failed': 'Failed', 'not_run': 'Not run'}


# Uzun bir toplamı okunur tut: kısa ise tek satır, değilse her terim ayrı satırda (aligned). Yalnız biçimdir; terimler aynıdır.
LINE_LIMIT = 110      # yaklaşık görünür karakter (LaTeX komutları sayılmaz); bunu aşan satır PDF'te okunmaz kadar küçülürdü


def visible_len(text):
    """LaTeX kaynağının yaklaşık görünür genişliği: komutlar, süslü parantezler ve boşluk komutları sayılmaz;
    \\frac{a}{b} iki katlı yazıldığı için payın ve paydanın uzun olanı kadar sayılır (kaba ama tutarlı ölçü)."""
    text = re.sub(r'\\frac\{([^{}]*)\}\{([^{}]*)\}', lambda m: max(m.group(1), m.group(2), key=len), text)
    text = re.sub(r'\\text\{([^{}]*)\}', r'\1', text)
    return len(re.sub(r'\\[a-zA-Z]+|[{}]|\\[,;!]|\s{2,}', '', text))


def stack(lhs, terms, tail=''):
    one = ' + '.join(terms).replace('+ -', '- ')
    if visible_len(lhs + one + tail) <= LINE_LIMIT:
        return f'{lhs} = {one}{tail}'
    rows = [f'{lhs} &= {terms[0]}'] + [('&\\quad ' + (t if t.startswith('-') else '+ ' + t)) for t in terms[1:]]
    tail = tail.lstrip(',').lstrip()                  # ayrı satıra inen ek (ör. x aralığı) virgülle başlamasın
    return '\\begin{aligned} ' + ' \\\\ '.join(rows) + (f' \\\\ &{tail}' if tail else '') + ' \\end{aligned}'


def chunk_rows(terms, width=110):
    """İşaretli terim listesini (ör. ['a', '- b', '+ c']) en fazla width karakterlik satırlara toplar."""
    rows, current = [], ''
    for t in terms:
        if current and len(current) + len(t) + 1 > width:
            rows.append(current)
            current = t
        else:
            current = f'{current} {t}' if current else t
    return rows + ([current] if current else [])


def signed_terms(expr):
    out = []
    args = sp.Add.make_args(sp.sympify(expr))
    try:                                              # azalan kuvvet sırası (tek satırlık LaTeX ile aynı düzen)
        args = sorted(args, key=lambda t: -sp.Poly(t, s).total_degree() if t.free_symbols <= {s} else 0)
    except Exception:
        pass
    for k, t in enumerate(args):
        body = sp.latex(t)
        out.append(body if k == 0 else (body if body.startswith('-') else '+ ' + body))
    return out


def sum_terms(expr):
    return [sp.latex(t) for t in sp.Add.make_args(expr)]


def product_parts(terms):
    """product_tex ile aynı çarpımlar, işaretleriyle ayrı ayrı (satıra bölmek için)."""
    out = []
    for sign, factors in terms:
        if any(g == 0 for g in factors):
            continue
        body = ''.join(f'\\left({sp.latex(g)}\\right)' for g in factors if g != 1) or '1'
        out.append(('- ' if sign < 0 else '+ ' if out else '') + body)
    return out or ['0']


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


# Motorun integral kaydını kâğıt satırına çevirir: yerine konmuş çarpım → açılmış integrand → ilkel fonksiyon → değer.
def format_integral(step):
    upper = sp.latex(step['h'])
    raw = f"{step['label']} &= \\int_{{0}}^{{{upper}}}\\left[{product_tex(step['terms'])}\\right]ds"
    wrap = lambda rest: f'\\begin{{aligned}} {raw} \\\\ &{rest} \\end{{aligned}}'
    if step['method'] == 'zero':
        return wrap('= 0')
    if step['method'] == 'polynomial':
        expanded = f"= \\int_{{0}}^{{{upper}}}\\left({sp.latex(step['integrand'])}\\right)ds"
        primitive = f"\\left[{sp.latex(step['primitive'])}\\right]_{{0}}^{{{upper}}} = {tex_value(step['value'])}"
        if visible_len(raw) <= LINE_LIMIT and visible_len(expanded + primitive) <= LINE_LIMIT:
            return wrap(f'{expanded} = {primitive}')
        # Uzun (yüksek dereceli) integral: her adım kendi satır(lar)ında; terimler aynen korunur, yalnız satıra bölünür.
        product = chunk_rows(product_parts(step['terms']))
        rows = [f"{step['label']} &= \\int_{{0}}^{{{upper}}}\\Big[ {product[0]}"] + [f'&\\qquad {r_}' for r_ in product[1:]]
        rows[-1] += ' \\Big]\\,ds'
        integrand = chunk_rows(signed_terms(step['integrand']))
        rows += [f'&= \\int_{{0}}^{{{upper}}}\\Big( {integrand[0]}'] + [f'&\\qquad {r_}' for r_ in integrand[1:]]
        rows[-1] += ' \\Big)\\,ds'
        prim = chunk_rows(signed_terms(step['primitive']))
        rows += [f'&= \\Big[ {prim[0]}'] + [f'&\\qquad {r_}' for r_ in prim[1:]]
        rows[-1] += f' \\Big]_{{0}}^{{{upper}}}'
        rows.append(f"&= {tex_value(step['value'])}")
        return '\\begin{aligned} ' + ' \\\\ '.join(rows) + ' \\end{aligned}'
    if step['method'] == 'symbolic':
        if step.get('primitive') is None:
            return wrap(f"= {tex_value(step['value'])}")
        expanded = f"= \\int_{{0}}^{{{upper}}}\\left({sp.latex(step['integrand'])}\\right)ds"
        primitive = f"\\left[{sp.latex(step['primitive'])}\\right]_{{0}}^{{{upper}}}"
        rows = chunk_rows(signed_terms(step['primitive']), 80)   # polinom olmayan ilkel: her adım ayrı satırda
        prim = '&= \\Big[ ' + ' \\\\ &\\qquad '.join(rows) + f' \\Big]_{{0}}^{{{upper}}}'
        return wrap(f"{expanded} \\\\ {prim} \\\\ &= {tex_value(step['value'])}")
    return wrap(f"\\approx {dec(step['value'])}")


# Motorun şekil fonksiyonu kurulum kaydını (Vandermonde yolu) blok blok yazar.
def format_shape(c):
    degree, h, cs, hats = c['degree'], c['h'], c['cs'], c['hats']
    B = []
    B.append(f'u^{{e}}(s) = {sp.latex(c["poly"])},\\qquad 0 \\le s \\le h = {sp.latex(h)}\\qquad '
             f'\\text{{(polynomial of degree {degree}: {degree + 1} coefficients, {degree + 1} nodes)}}')
    B.append('s_k = (k - 1)\\,\\frac{h}{p}:\\qquad ' + ',\\ '.join(f's_{{{k + 1}}} = {sp.latex(v)}' for k, v in enumerate(c['nodes'])))
    B.append('\\begin{aligned} ' + ' \\\\ '.join(f'u^{{e}}(s_{{{k + 1}}}) &= {sp.latex(cond)} = {sp.latex(hats[k])}'
                                               for k, cond in enumerate(c['conds'])) + ' \\end{aligned}\\qquad \\text{(interpolation conditions)}')
    B.append(f'{tex_matrix(c["V"].tolist())}\\,{tex_column(cs)} = {tex_column(hats)}\\qquad \\mathbf{{V}}\\mathbf{{c}} = \\mathbf{{a}}^{{(e)}}'
             f'\\quad \\text{{(interpolation matrix }} \\mathbf{{V}}\\text{{, }} V_{{km}} = s_k^{{\\,m}}\\text{{)}}')
    B.append(f's_1 = 0 \\;\\Rightarrow\\; c_0 = {sp.latex(c["sol"][cs[0]])}')
    B.append(f'{tex_matrix(c["W"].tolist())}\\,{tex_column(cs[1:])} = {tex_column(c["rhs"])}\\qquad \\text{{(remaining {degree} equations, }} c_0 \\text{{ substituted)}}')
    B.append('\\begin{aligned} ' + ' \\\\ '.join(f'{sp.latex(cm)} &= {sp.latex(sp.collect(sp.expand(c["sol"][cm]), hats))}' for cm in cs)
             + ' \\end{aligned}\\qquad \\text{(solution)}')
    B.append(stack('u^{e}(s)', [f'{sp.latex(hv)}\\left({sp.latex(Nk)}\\right)' for hv, Nk in zip(hats, c['shapes'])],
                   f' = \\sum_{{k=1}}^{{{degree + 1}}} a^{{(e)}}_k\\,N_k(s)\\qquad \\text{{(grouped by the coefficients }}a^{{(e)}}_k\\text{{)}}'))
    B.append('\\begin{aligned} ' + ' \\\\ '.join(f"N_{{{k + 1}}}(s) &= {sp.latex(Nk)}, &\\qquad N_{{{k + 1}}}'(s) &= {sp.latex(dNk)}"
                                               for k, (Nk, dNk) in enumerate(zip(c['shapes'], c['derivatives']))) + ' \\end{aligned}')
    B.append(f'\\big[N_i(s_j)\\big] = {tex_matrix(c["delta"].tolist())} = \\mathbf{{I}}\\quad (N_i(s_j) = \\delta_{{ij}}),\\qquad '
             f'\\sum_k N_k(s) = {sp.latex(c["partition"])}\\ \\text{{(partition of unity)}}')
    B.append(f'\\sum_k N_k(s)\\,s_k^{{\\,m}} = s^m,\\quad m = 0, \\dots, {degree}\\qquad \\text{{(polynomial completeness: every polynomial up to degree {degree} is represented exactly)}}')
    return B


# Testlerin ve dışarıdan çağıranların kolaylığı için: kurulumu motordan alıp biçimler (kâğıt kendisi r['shape']'i kullanır).
def shape_function_derivation(degree, h, N_engine):
    return format_shape(shape_function_construction(degree, h, N_engine))


# Sürekli sınıflandırma (Katman 4) ve katsayı dejenereliği (Katman 2) kanıt düzeyiyle yazılır.
def classification_lines(r, c):
    lines, g = [], r.get('degeneracy')
    lead = '\\alpha' if r['order'] == 2 else '\\beta'
    if g and g['kind'] != 'uniform':
        where = ',\\ '.join(sp.latex(pt) if g['method'] == 'symbolic' else dec(pt, 6) for pt in g['points'])
        name = 'sign-changing' if g['kind'] == 'sign_changing' else 'dejenere'
        lines.append(f'{lead}(x^*) = 0,\\ x^* \\in \\{{{where}\\}}\\qquad \\text{{({name} coefficient; {EVIDENCE_LABEL[g["method"]]} detection; no uniform ellipticity, standard convergence theory not guaranteed)}}')
    if c:
        if c['kind'] == 'unique':
            extra = (';\\ ' + '; '.join(c['assumptions'])) if c['assumptions'] else ''
            lines.append(f'\\text{{Continuous problem: unique solution ({EVIDENCE_LABEL.get(c["method"], c["method"])}; {c["evidence"]}{extra})}}')
        elif c['kind'] in ('no_solution', 'non_unique'):          # yalnız sayısal kanıtlı olanlar buraya ulaşır
            name = 'no solution' if c['kind'] == 'no_solution' else 'not unique'
            lines.append(f'\\text{{Continuous problem: {name} ({EVIDENCE_LABEL.get(c["method"], c["method"])}; not a proof; {c["evidence"]})}}')
        else:
            lines.append(f'\\text{{Continuous problem not classified: {c["reason"]}}}')
    return lines


def wrap_text(text, width):
    """Uzun düz metni ', ' ve '; ' sınırlarından width karakteri aşmayan satırlara böler (kâğıtta okunur kalsın)."""
    lines, current = [], ''
    for piece in re.split(r'(?<=[,;])\s+', text):
        if current and len(current) + 1 + len(piece) > width:
            lines.append(current)
            current = piece
        else:
            current = f'{current} {piece}' if current else piece
    return lines + [current] if current else lines


def text_block(text, width=95):
    chunks = wrap_text(text, width)
    if len(chunks) == 1:
        return f'\\text{{{chunks[0]}}}'
    return '\\begin{aligned} ' + ' \\\\ '.join(f'&\\text{{{c}}}' for c in chunks) + ' \\end{aligned}'


def raw_tex(text):
    """Kullanıcının yazdığı metin, karakteri karakterine (\\texttt). LaTeX'e özel karakterler kaçırılır; ^ build_tex'te çevrilir."""
    if not text:
        return '\\text{(empty)}'
    return '\\texttt{' + ''.join('\\' + ch if ch in '{}%&#$_' else ch for ch in text) + '}'


def build_sheet(valid, r, V):
    """Kâğıt: ders kitabındaki matematiksel sırayla ve yalnız hesaplanmış veriden. Verilenler (girdi → iç gösterim izi,
    güçlü form, koşullar, analitik çözüm durumu) → gösterim → ağ → şekil fonksiyonları → yaklaşık çözüm → Galerkin ağırlıklı
    kalıntı → (yalnız 2. mertebe) zayıf form → eleman integralleri → montaj → sınır koşulları → çözüm → yeniden kurulum →
    hata ölçüleri → doğrulama ve tanılar → uyarılar ve sınırlar."""
    alpha, beta, gamma, f = r['alpha'], r['beta'], r['gamma'], r['f']
    x0, x1, n, degree, h, dofs, order = r['x0'], r['x1'], r['n'], r['degree'], r['h'], r['dofs'], r['order']
    form = r['formulation']
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
            return f'u({P}) = {sp.latex(value)}\\ \\text{{(essential condition)}}'
        if kind == 'neumann':
            return f"u'({P}) = {sp.latex(value)}\\ \\text{{(derivative condition; flux }} Q = \\alpha({P})\\,u'({P})\\text{{)}}"
        if kind == 'flux':
            return f"\\alpha({P})\\,u'({P}) = {sp.latex(value)}\\ \\text{{(flux condition; }} Q = {sp.latex(value)}\\text{{)}}"
        if kind == 'both':
            return f"u({P}) = {sp.latex(value)},\\ u'({P}) = {sp.latex(r['slopes'][side])}\\ \\text{{(initial conditions)}}"
        if r['ivp'] or order != 2:
            return None
        if alpha.subs(x, point) == 0:
            return f"\\alpha({P}) = 0\\ \\text{{(not given: the flux term vanishes by itself; not an independent condition)}}"
        return f"\\alpha u'\\big|_{{x={P}}} = 0\\ \\text{{(not given: natural condition)}}"

    exact = V['exact']
    solution = exact['u']
    level = {'verified_symbolic': 'verified symbolically', 'verified_numeric': 'verified numerically (50 digits; not a proof)'}

    # --- Verilenler
    lhs_tex = tex_equation(alpha, beta, gamma, f).split(' = ')[0]
    add('## Given')
    mark = {True: '\\checkmark', False: '\\times'}
    rows = ' \\\\ '.join(f"\\text{{{t['item']}}} & {raw_tex(t['raw'])} & {t['internal']} & {mark[t['ok']]}\\ \\text{{{EVIDENCE_LABEL.get(t['how'], t['how'])}}}"
                         for t in valid['trace'])
    D(f'\\begin{{array}}{{l|l|l|l}} \\text{{field}} & \\text{{input (verbatim)}} & \\text{{internal form}} & \\text{{equivalence}} \\\\ \\hline {rows} \\end{{array}}')
    D(tex_equation(alpha, beta, gamma, f) + f', \\qquad x \\in \\Omega = ({X0},\\, {X1})')
    D(f'\\mathcal{{L}}(u) = {lhs_tex}\\qquad \\text{{(differential operator of order {order}; strong form: }} \\mathcal{{L}}(u) = f\\text{{)}}')
    D(',\\qquad '.join(c for c in (condition_tex('left', x0), condition_tex('right', x1)) if c))
    D(f'n = {n},\\qquad h = \\frac{{x_1 - x_0}}{{n}} = {sp.latex(h)},\\qquad p = {degree}\\ \\text{{({DEGREE_NAMES[degree]} element: {degree + 1} nodes per element)}}')
    D(f'N = np + 1 = {dofs}\\ \\text{{nodes (nodal coefficient }} a_j\\text{{); unknowns after essential conditions: }} {len(r["free"])}')
    method_text = ('Bubnov–Galerkin' if not r['ivp'] else 'Petrov–Galerkin (test functions vanish at the end with unknown flux)')
    D(f'\\text{{Method: {method_text}}}')
    for line in classification_lines(r, V['classification']):
        D(line)
    if solution is not None:
        D(f'u(x) = {sp.latex(solution)}\\qquad \\text{{(analytical solution; {level[exact["status"]]}: equation and all conditions)}}')
    elif exact['status'] == 'rejected':
        failed = next(c for c in exact['checks'] if not c[1])
        D(f"\\text{{analytical candidate }} {sp.latex(exact['candidate'])}\\ \\text{{rejected: }} {failed[0][1]}\\ \\text{{is not satisfied (not used)}}")
    else:
        D(f'\\text{{no analytical solution: {exact["reason"]}}}')

    section('Notation')
    D('\\begin{aligned} u(x)&:\\ \\text{exact solution (if known)} & u_h(x)&:\\ \\text{finite element approximation (global)} \\\\ '
      'N_j(x)&:\\ \\text{global shape function (node } j\\text{)} & u_h^{(e)}(x)&:\\ \\text{local approximation on element } e \\\\ '
      'N_k(s),\\ N_k^{(e)}(x)&:\\ \\text{local shape (interpolation) function} & \\phi_i(x)&:\\ \\text{Galerkin weighting function} \\\\ '
      'a_j&:\\ \\text{unknown global nodal coefficient} & a^{(e)}_k&:\\ \\text{local coefficient of element } e \\\\ '
      'c_m&:\\ \\text{coefficient of the element polynomial} & R(u_h, x)&:\\ \\text{residual} \\\\ '
      'K_{ij},\\ F_i&:\\ \\text{system matrix and load vector entries} & Q&:\\ \\text{end flux } \\alpha u\' \\end{aligned}'
      '\\qquad \\text{(Heinrich \\& Pepper notation)}')
    D('N_j(x) = \\begin{cases} N_k^{(e)}(x), & x \\in \\Omega_e,\\ j = (e-1)p + k \\\\ 0, & \\text{on the other elements} \\end{cases},\\qquad '
      'a^{(e)}_k = a_{(e-1)p+k},\\qquad \\phi_i = N_i\\ \\text{(Bubnov–Galerkin)}')

    section('Mesh')
    nodes = form['nodes']
    D('\\big[x_j\\big] = ' + tex_row(nodes) + f'\\qquad \\Omega_e = [x_e,\\ x_e + h],\\ e = 1, \\dots, {n}')
    shown_e = list(range(n)) if n <= 6 else [0, 1, None, n - 1]
    table = ' \\\\ '.join('\\vdots & \\vdots & \\vdots' if e is None else
                          f'{e + 1} & [{sp.latex(x0 + e * h)},\\ {sp.latex(x0 + (e + 1) * h)}] & '
                          + ',\\ '.join(f'{e * degree + k + 1}' for k in range(degree + 1)) for e in shown_e)
    D(f'\\begin{{array}}{{c|c|c}} e & \\Omega_e & \\text{{local }} k = 1, \\dots, {degree + 1} \\;\\mapsto\\; \\text{{global }} j = (e-1)p + k \\\\ \\hline {table} \\end{{array}}')

    section('Shape functions (from polynomial interpolation)')
    D('\\text{Reference element: } s \\in [0, h],\\quad s = x - x_e\\qquad \\text{(every } \\Omega_e \\text{ is this element shifted by } x = x_e + s\\text{; Jacobian } dx/ds = 1\\text{; the derivation is the same on every element)}')
    for block in format_shape(r['shape']):
        D(block)
    D('x = x_e + s,\\quad dx = ds,\\quad \\frac{dN_k^{(e)}}{dx} = \\frac{dN_k}{ds}\\ (\\text{Jacobian } 1),\\qquad N_k^{(e)}(x) = N_k(x - x_e),'
      '\\qquad a^{(e)}_k = a_{(e-1)p+k}\\qquad \\text{(local → global)}')
    D('\\text{continuity: neighbouring elements share the coefficient } a_j \\text{ at the common node} \\Rightarrow u_h \\in C^0;\\ '
      'u_h\' \\text{ is in general discontinuous at the nodes}')
    for e in form['shown']:
        E = form['elements'][e]
        shapes = [f'N_{{{k + 1}}}^{{({e + 1})}}(x) = {sp.latex(Nk)}' for k, Nk in enumerate(E['N'])]
        mapping = ',\\ '.join(f'a^{{({e + 1})}}_{{{k + 1}}} = {sp.latex(a)}' for k, a in enumerate(E['coefficients']))
        if visible_len(''.join(shapes) + mapping) <= LINE_LIMIT:
            D(f'e = {e + 1}:\\qquad ' + ',\\quad '.join(shapes) + ',\\qquad ' + mapping)
        else:
            D(f'e = {e + 1}:\\qquad \\begin{{aligned}} ' + ' \\\\ '.join(t.replace(' = ', ' &= ', 1) for t in shapes) + ' \\end{aligned}\\qquad ' + mapping)

    section('Approximate solution')
    D(f'u_h(x) = \\sum_{{j=1}}^{{{dofs}}} a_j\\,N_j(x)\\qquad \\text{{(}} a_j \\text{{: nodal coefficients to be determined)}}')
    D('N_j(x_i) = \\delta_{ij} \\;\\Rightarrow\\; u_h(x_i) = a_i\\qquad \\text{(the coefficients are the nodal values)}')
    D('x \\in \\Omega_e:\\quad u_h(x) = u_h^{(e)}(x) = \\sum_{k=1}^{p+1} a^{(e)}_k\\,N_k^{(e)}(x)')
    for e in form['shown']:
        E = form['elements'][e]
        D(stack(f'u_h^{{({e + 1})}}(x)', [f'{sp.latex(a)}\\left({sp.latex(Nk)}\\right)' for a, Nk in zip(E['coefficients'], E['N'])],
                f',\\qquad x \\in [{sp.latex(E["xe"])}, {sp.latex(E["xe"] + h)}]'))
    if n > len(form['shown']):
        D('\\text{the other elements have the same structure}')

    # --- Galerkin ağırlıklı kalıntı
    section('Galerkin Weighted-Residual Formulation')
    D(f'R(u_h, x) = \\mathcal{{L}}(u_h) - f = \\left[{lhs_tex}\\right]_{{u = u_h}}' + ('' if f == 0 else f' - \\left({sp.latex(f)}\\right)'))
    interior_name = 'R_e(x)' if order == 2 else 'R(u_h, x)\\big|_{\\Omega_{e}}'
    for e in form['shown']:
        label = f'R_{{{e + 1}}}(x)' if order == 2 else f'R(u_h, x)\\big|_{{\\Omega_{{{e + 1}}}}}'
        res, note = form['elements'][e]['residual'], ('\\qquad \\text{(classical residual inside the element)}' if order == 2 else '')
        D(f'{label} = {sp.latex(res)}{note}' if visible_len(sp.latex(res)) <= LINE_LIMIT else stack(label, sum_terms(res), note))
    if order == 2:
        D("\\mathcal{L}(u_h) - f = \\sum_e R_e\\,\\chi_{\\Omega_e} + \\sum_{k} J_k\\,\\delta(x - x_k),\\qquad "
          "J_k = \\alpha(x_k)\\big[u_h'(x_k^+) - u_h'(x_k^-)\\big]\\qquad \\text{(distributional residual; derivative jump at interior nodes)}")
        D('\\text{the element residual } R_e \\text{ alone is not the global residual; a } C^0 \\text{ element also has the nodal terms } J_k')
    D(f'\\int_{{{X0}}}^{{{X1}}} \\phi_i(x)\\,R(u_h, x)\\,dx = 0\\qquad \\text{{(the residual vanishes on average with respect to the weighting functions)}}')
    D('\\phi_i(x) = N_i(x)\\qquad \\text{(' + ('Bubnov–Galerkin: the weighting functions are the shape functions' if not r['ivp']
      else 'Petrov–Galerkin: the weights are chosen from the shape functions, but the test space differs from the trial space') + ')}')
    tests_used = list(r['rows'])
    if tests_used:
        D(f'\\phi_i \\in \\{{N_i\\}},\\quad {index_range(tests_used[0] + 1, tests_used[-1] + 1)}\\qquad \\text{{({len(tests_used)} equations, {len(r["free"])} unknowns)}}')
    else:                                     # tüm düğümler temel koşulla belirli: ağırlık uzayı boş, çözülecek denklem yok
        D('\\text{empty weighting space: all nodal values are fixed by essential conditions (0 equations, 0 unknowns)}')
    for i in range(dofs):
        if i in tests_used:
            continue
        P = 'x_0' if i == 0 else 'x_1'
        if r['ivp']:
            D(f'\\phi_{{{i + 1}}} = N_{{{i + 1}}}\\ \\text{{removed: the flux at }} {P}\\text{{ is unknown, the weighting function is chosen to vanish there}}')
        else:
            D(f'\\phi_{{{i + 1}}} = N_{{{i + 1}}}\\ \\text{{removed: }} u({P})\\ \\text{{is given (essential condition) }}\\Rightarrow a_{{{i + 1}}}\\ \\text{{known; weighting space: }} \\phi({P}) = 0')

    def matrix_definitions():
        D(f'\\mathbf{{a}} = \\begin{{bmatrix}} a_1 & \\cdots & a_{{{dofs}}} \\end{{bmatrix}}^{{\\mathsf{{T}}}} \\in \\mathbb{{R}}^{{{dofs}}},\\quad '
          f'\\mathbf{{K}} \\in \\mathbb{{R}}^{{{dofs}\\times{dofs}}},\\quad \\mathbf{{F}} \\in \\mathbb{{R}}^{{{dofs}}}\\qquad '
          '\\text{(row } i \\leftrightarrow \\phi_i\\text{, column } j \\leftrightarrow a_j\\text{; } K_{ij}\\text{: coefficient of } a_j\\text{)}')
        D('\\sum_{j} K_{ij}\\,a_j = F_i \\;\\Longleftrightarrow\\; \\mathbf{K}\\,\\mathbf{a} = \\mathbf{F}')

    if order == 1:
        D('\\text{First order: } R \\text{ contains only } u_h\'\\text{, which is integrable for a } C^0 \\text{ approximation} \\Rightarrow \\text{no integration by parts is needed}')
        for row in form['galerkin_rows']:
            parts = ' + '.join(f'\\int_{{{sp.latex(p_["xe"])}}}^{{{sp.latex(p_["xe"] + h)}}}\\left({sp.latex(p_["weight"])}\\right)\\left[{sp.latex(p_["residual"])}\\right]dx'
                               for p_ in row['parts'])
            if visible_len(parts + sp.latex(row['total'])) <= LINE_LIMIT:
                D(f'\\phi_{{{row["i"] + 1}}}:\\quad {parts} = {sp.latex(row["total"])} = 0')
            else:                                         # uzun: ağırlık ve kalıntı yukarıda yazıldı, adlarıyla anılır
                k_of = lambda p_: row['i'] - p_['e'] * degree + 1
                named = ' + '.join(f'\\int_{{{sp.latex(p_["xe"])}}}^{{{sp.latex(p_["xe"] + h)}}} N_{{{k_of(p_)}}}^{{({p_["e"] + 1})}}(x)\\,'
                                   f'R(u_h, x)\\big|_{{\\Omega_{{{p_["e"] + 1}}}}}\\,dx' for p_ in row['parts'])
                total = chunk_rows(signed_terms(row['total']))
                D(f'\\phi_{{{row["i"] + 1}}}:\\quad \\begin{{aligned}} &{named} \\\\ &= ' + ' \\\\ &\\quad '.join(total) + ' = 0 \\end{aligned}')
        k_int = tex_terms([(beta, "\\phi_i\\,N_j'"), (gamma, '\\phi_i\\,N_j')])
        f_int = '0' if f == 0 else f'\\int_{{{X0}}}^{{{X1}}} \\phi_i\\left({sp.latex(f)}\\right)dx'
        D(f'u_h = \\sum_j a_j N_j:\\quad \\int_{{{X0}}}^{{{X1}}} \\phi_i\\,R\\,dx = \\sum_j a_j \\underbrace{{\\int_{{{X0}}}^{{{X1}}} \\left({k_int}\\right)dx}}_{{K_{{ij}}}}'
          f' - \\underbrace{{{f_int}}}_{{F_i}} = 0,\\qquad \\phi_i = N_i')
        matrix_definitions()
    else:
        D("\\text{Second order: } R \\ni (\\alpha u_h')' \\text{; for a } C^0 \\text{ approximation } u_h' \\text{ is discontinuous at the nodes, so } (\\alpha u_h')' \\text{ is undefined there in the classical sense}")
        D("\\Rightarrow \\text{one derivative is moved onto the weighting function (next section)}")
        D(f"\\int_{{{X0}}}^{{{X1}}} \\phi_i\\left[{lhs_tex.replace('u', 'u_h')}" + ('' if f == 0 else f' - \\left({sp.latex(f)}\\right)') + '\\right]dx = 0')

        section('Weak Formulation (Integration by Parts)')
        D(f"\\int_{{{X0}}}^{{{X1}}} \\phi_i\\,(\\alpha u_h')'\\,dx = \\Big[\\phi_i\\,\\alpha u_h'\\Big]_{{{X0}}}^{{{X1}}} - \\int_{{{X0}}}^{{{X1}}} \\phi_i'\\,\\alpha u_h'\\,dx"
          "\\qquad \\text{(integration by parts)}")
        D("\\text{Result: only the first derivative of } u_h \\text{ is required (derivative order 2 → 1); a } C^0 \\text{ approximation suffices; the end flux } Q = \\alpha u' \\text{ appears}")
        weak_left = tex_terms([(alpha, "\\phi_i'u_h'"), (-beta, "\\phi_i\\,u_h'"), (-gamma, '\\phi_i\\,u_h')])
        load = '' if f == 0 else f" - \\int_{{{X0}}}^{{{X1}}} \\phi_i\\left({sp.latex(f)}\\right)dx"
        D(f"\\int_{{{X0}}}^{{{X1}}}\\left({weak_left}\\right)dx = \\phi_i({X1})\\,Q_1 - \\phi_i({X0})\\,Q_0{load},\\qquad Q = \\alpha u'")
        for side, point, name in (('left', x0, 'Q_0'), ('right', x1, 'Q_1')):
            kind, value = r['conditions'][side]
            P = sp.latex(point)
            if kind == 'dirichlet':
                why = f'u({P})\\ \\text{{is given (essential condition): in the weighting space }} \\phi_i({P}) = 0 \\Rightarrow \\text{{the term is absent from the equations used; }} {name} \\text{{ is unknown and is computed after the solve}}'
            elif kind == 'neumann':
                why = f"{name} = \\alpha({P})\\,u'({P}) = {tex_value(r['Q'][side])}\\ \\text{{(derivative condition, natural: enters the flux multiplied by α; into the load vector)}}"
            elif kind == 'flux':
                why = f"{name} = \\alpha({P})\\,u'({P}) = {tex_value(r['Q'][side])}\\ \\text{{(flux condition, natural: enters the load vector directly)}}"
            elif kind == 'both':
                why = f"{name} = \\alpha({P})\\,u'({P}) = {tex_value(r['Q'][side])}\\ \\text{{(initial condition: the flux is known)}}"
            elif r['ivp']:
                why = f'{name}\\ \\text{{unknown: in the weighting space }} \\phi_i({P}) = 0 \\text{{ is chosen (Petrov–Galerkin), the term drops}}'
            elif alpha.subs(x, point) == 0:
                why = f'\\alpha({P}) = 0 \\Rightarrow {name} = 0\\ \\text{{(degenerate end)}}'
            else:
                why = f'{name} = 0\\ \\text{{(no condition given: natural condition }} \\alpha u\' = 0\\text{{)}}'
            D(f'x = {P}:\\quad ' + why)
        k_int = tex_terms([(alpha, "\\phi_i'N_j'"), (-beta, "\\phi_i\\,N_j'"), (-gamma, '\\phi_i\\,N_j')])
        f_int = '' if f == 0 else f' - \\int_{{{X0}}}^{{{X1}}} \\phi_i\\left({sp.latex(f)}\\right)dx'
        D(f'u_h = \\sum_j a_j N_j,\\ \\phi_i = N_i:\\quad \\sum_j a_j \\underbrace{{\\int_{{{X0}}}^{{{X1}}} \\left({k_int}\\right)dx}}_{{K_{{ij}}}}'
          f' = \\underbrace{{\\phi_i({X1})\\,Q_1 - \\phi_i({X0})\\,Q_0{f_int}}}_{{F_i = F_{{b,i}} + F_{{l,i}}}}')
        matrix_definitions()

    # --- eleman integralleri
    section('Element matrices and load vectors')
    D(f'K_{{ij}} = \\sum_e K^{{(e)}}_{{ij}},\\qquad F_{{l,i}} = \\sum_e F^{{(e)}}_{{l,i}}\\qquad \\text{{(}}\\phi_i, N_j \\text{{ are nonzero only on the elements containing their nodes)}}')
    D(f'i, j = 1, \\dots, {degree + 1}\\ \\text{{(local)}},\\qquad \\mathbf{{K}}^{{(e)}} \\in \\mathbb{{R}}^{{{degree + 1}\\times{degree + 1}}},\\ \\mathbf{{F}}^{{(e)}} \\in \\mathbb{{R}}^{{{degree + 1}}}')
    integrand_x = (tex_terms([(alpha, "N_i^{(e)\\prime}N_j^{(e)\\prime}"), (-beta, "N_i^{(e)}N_j^{(e)\\prime}"), (-gamma, 'N_i^{(e)}N_j^{(e)}')]) if order == 2
                   else tex_terms([(beta, "N_i^{(e)}N_j^{(e)\\prime}"), (gamma, 'N_i^{(e)}N_j^{(e)}')]))
    at_s = lambda g: g.subs(x, sp.Symbol('x_e') + s)          # katsayı yerel koordinatta: x = x_e + s
    integrand_K = (tex_terms([(at_s(alpha), "N_i'N_j'"), (-at_s(beta), "N_i\\,N_j'"), (-at_s(gamma), 'N_i\\,N_j')]) if order == 2
                   else tex_terms([(at_s(beta), "N_i\\,N_j'"), (at_s(gamma), 'N_i\\,N_j')]))
    D(f'K^{{(e)}}_{{ij}} = \\int_{{x_e}}^{{x_e + h}}\\left({integrand_x}\\right)dx,\\qquad '
      f'F^{{(e)}}_{{l,i}} = {"-" if order == 2 else ""}\\int_{{x_e}}^{{x_e + h}} N_i^{{(e)}}(x)\\,f(x)\\,dx'
      + ('' if order == 2 else '\\qquad (\\mathbf{F} = \\mathbf{F}_l\\text{; no boundary term})'))
    D('x = x_e + s,\\quad J = \\frac{dx}{ds} = 1,\\quad dx = J\\,ds = ds,\\quad \\frac{dN_i^{(e)}}{dx} = \\frac{1}{J}\\frac{dN_i}{ds} = N_i\'(s),'
      '\\quad x \\in [x_e, x_e + h] \\mapsto s \\in [0, h]\\qquad \\text{(coordinate transformation)}')
    D(f'K^{{(e)}}_{{ij}} = \\int_0^h\\left({integrand_K}\\right)ds,\\qquad F^{{(e)}}_{{l,i}} = {"-" if order == 2 else ""}\\int_0^h N_i(s)\\,f(x_e + s)\\,ds'
      + "\\qquad (N_i' = dN_i/ds)")
    methods = {step['method'] for record in r['elements'] for step in record['K_steps'] + record['F_steps']}
    if not r['symbolic']:
        D('\\text{Symbolic time limit exceeded: all integrals numerical (mpmath.quad, adaptive tanh–sinh rule, 15 digits); steps not written}')
    elif 'numeric' in methods:
        D('\\approx:\\ \\text{no closed form found; numerical integral (mpmath.quad, adaptive tanh–sinh rule, 15 digits)}')
    else:
        D('\\text{The integrals shown were computed exactly (symbolically) with SymPy.}')
    if n > DETAIL_LIMIT and r['symbolic']:
        D(f'\\text{{elements not shown (e = 3, …, {n - 1}): numerical integral (mpmath.quad, tanh–sinh); values in the assembly}}')
    if r['uniform_K'] and n > 1:
        D('\\alpha, \\beta, \\gamma \\text{ constant},\\ h \\text{ equal on every element},\\ J = 1 \\Rightarrow \\text{the integrand in } s \\text{ does not depend on } e \\Rightarrow \\mathbf{K}^{(e)} = \\mathbf{K}^{(1)}\\ \\text{(all } e\\text{)}')
    if r['uniform_F'] and n > 1:
        D('f \\text{ constant},\\ h \\text{ equal},\\ J = 1 \\Rightarrow \\mathbf{F}^{(e)}_l = \\mathbf{F}^{(1)}_l\\ \\text{(all } e\\text{)}')
    for record in r['elements']:
        e = record['e']
        if e not in r['detailed']:
            continue
        has_K, has_F = not (r['uniform_K'] and e > 0), not (r['uniform_F'] and e > 0)
        if not (has_K or has_F):
            continue
        xe = record['xe']
        D(f'e = {e + 1}:\\qquad x_e = {sp.latex(xe)},\\qquad \\int_{{{sp.latex(xe)}}}^{{{sp.latex(xe + h)}}}(\\cdots)\\,dx = '
          f'\\int_0^{{{sp.latex(h)}}}(\\cdots)\\big|_{{x = {sp.latex(xe + s)}}}\\,ds')
        if not has_K:
            D(f'\\mathbf{{K}}^{{({e + 1})}} = \\mathbf{{K}}^{{(1)}} = {tex_matrix(record["K"].tolist())}\\qquad \\text{{(reason above)}}')
        if has_K:
            if r['symbolic']:
                for step in record['K_steps']:
                    D(format_integral(step))
            D(f'\\mathbf{{K}}^{{({e + 1})}} = {tex_matrix(record["K"].tolist())}')
        if has_F:
            if r['symbolic']:
                for step in record['F_steps']:
                    D(format_integral(step))
            D(f'\\mathbf{{F}}_l^{{({e + 1})}} = {tex_column(list(record["F"]))}')
    same = [record['e'] + 1 for record in r['elements'] if record['e'] in r['detailed'] and record['e'] > 0 and r['uniform_K'] and r['uniform_F']]
    if same:
        D(f'e = {index_range(same[0], same[-1])[4:]}:\\qquad \\mathbf{{K}}^{{(e)}} = \\mathbf{{K}}^{{(1)}},\\quad \\mathbf{{F}}_l^{{(e)}} = \\mathbf{{F}}_l^{{(1)}}\\qquad \\text{{(same integrand, same limits; reason above)}}')
    if order == 2 and not alpha.free_symbols and alpha < 0:
        first_K, first_F = r['elements'][0]['K'], r['elements'][0]['F']
        D(f'\\alpha = {sp.latex(alpha)} < 0:\\qquad \\tilde{{\\mathbf{{K}}}}^{{(1)}} = -\\mathbf{{K}}^{{(1)}} = {tex_matrix((-first_K).tolist())},'
          f'\\qquad \\tilde{{\\mathbf{{F}}}}_l^{{(1)}} = -\\mathbf{{F}}_l^{{(1)}} = {tex_column(list(-first_F))}')
        D('\\tilde{\\mathbf{K}}\\,\\mathbf{a} = \\tilde{\\mathbf{F}} \\iff (-\\mathbf{K})\\,\\mathbf{a} = -\\mathbf{F} \\iff '
          '\\mathbf{K}\\,\\mathbf{a} = \\mathbf{F}\\qquad \\text{(same }\\mathbf{a}\\text{; with the standard positive stiffness sign)}')

    # --- montaj
    section('Assembly')
    D('K_{IJ} = \\sum_{e} K^{(e)}_{ij},\\qquad F_I = \\sum_e F^{(e)}_i,\\qquad I = (e-1)p + i,\\ J = (e-1)p + j\\qquad \\text{(local } i,j \\to \\text{ global } I,J\\text{)}')
    for c in form['contributions']:
        terms = ' + '.join(f'K^{{({e + 1})}}_{{{i + 1}{j + 1}}}' for e, i, j, _ in c['parts'])
        values = ' + '.join(f'\\left({tex(v)}\\right)' for *_, v in c['parts'])
        D(f'K_{{{c["I"] + 1}{c["J"] + 1}}} = {terms} = {values} = {tex(c["value"])}' + (f'\\qquad \\text{{(shared node {c["I"] + 1})}}' if c['I'] == c['J'] else ''))
    for c in form['load_contributions']:
        D(f'F_{{l,{c["I"] + 1}}} = ' + ' + '.join(f'F^{{({e + 1})}}_{{l,{i + 1}}}' for e, i, _ in c['parts']) + ' = '
          + ' + '.join(f'\\left({tex(v)}\\right)' for *_, v in c['parts']) + f' = {tex(c["value"])}')
    if form['contributions']:
        D('\\text{the other entries come from a single element: } K_{IJ} = K^{(e)}_{ij},\\ \\text{empty entries: nodes without a common element}')
    D('\\mathbf{K} = ' + tex_matrix(r['K'].tolist()) + ',\\qquad \\mathbf{F}_l = ' + tex_column(list(r['F_l'])))
    if order == 2:
        D("\\mathbf{F}_b = \\begin{bmatrix} -Q_0 & 0 & \\cdots & 0 & Q_1 \\end{bmatrix}^{\\mathsf{T}}\\qquad \\text{(end nodes only; at interior nodes } \\phi_i(x_0) = \\phi_i(x_1) = 0\\text{)}")
        for row in form['row_equations']:
            D(f'\\phi_{{{row["i"] + 1}}}:\\quad {sp.latex(row["lhs"])} = {sp.latex(row["rhs"])}')

    # --- sınır koşulları
    section('Boundary conditions and reduced system')
    known = []
    for side, index, name in (('left', 0, 'Q_0'), ('right', dofs - 1, 'Q_1')):
        if index in r['fixed']:
            known.append(f'a_{{{index + 1}}} = {sp.latex(r["fixed"][index])}\\ \\text{{(essential)}}')
        if order == 2 and r['Q'][side] is not None:
            origin = 'from the initial condition' if r['conditions'][side][0] == 'both' else 'natural'
            known.append(f'{name} = {tex_value(r["Q"][side])}\\ \\text{{({origin})}}')
    if known:
        D(',\\qquad '.join(known))
    removed_rows = [i + 1 for i in range(dofs) if i not in r['rows']]
    removed_cols = [i + 1 for i in r['fixed']]
    if removed_rows:
        D('\\text{removed rows: } ' + ', '.join(str(i) for i in removed_rows) + '\\ (\\phi_i \\text{ not in the weighting space})')
    if removed_cols:
        D('\\text{removed columns: } ' + ', '.join(str(j) for j in removed_cols) + '\\ (a_j \\text{ known} \\Rightarrow K_{ij}a_j \\text{ to the right-hand side)}')
        D('F_{\\mathrm{red},i} = F_i - \\sum_{j\\ \\text{bilinen}} K_{ij}\\,a_j')
    if r['free']:
        unknowns = [sp.Symbol(f'a_{{{i + 1}}}') for i in r['free']]
        D(f'\\mathbf{{K}}_{{\\mathrm{{red}}}}\\,\\mathbf{{a}}_{{\\mathrm{{unknown}}}} = \\mathbf{{F}}_{{\\mathrm{{red}}}}:\\qquad {tex_matrix(r["K_ff"].tolist())}\\,{tex_column(unknowns)} = {tex_column(list(r["F_f"]))}')

    # --- çözüm
    section('Solution')
    D('\\mathbf{a} = ' + tex_row(r['a']) + (f"\\qquad \\text{{(decimal solution: {r['solve_note']})}}" if r.get('solve_note') else ''))
    if r['exact'] and any(not sp.sympify(v).is_Rational for v in r['a']):
        D('\\mathbf{a} \\approx ' + tex_row(r['a'], lambda v: dec(sp.sympify(v).evalf())))
    for side, name, point in (('left', 'Q_0', x0), ('right', 'Q_1', x1)):
        if side in r['reactions']:
            flux = r['reactions'][side]
            slope = (f",\\qquad u_h'({sp.latex(point)}) = \\frac{{{name}}}{{\\alpha({sp.latex(point)})}} = {tex_value(r['reaction_slopes'][side])}"
                     if side in r['reaction_slopes'] else '')
            D(f'{name} = {tex_value(flux)}{slope}\\qquad \\text{{(from the equation of the removed row)}}')
    if r['first_order']:
        fo = r['first_order']
        rd, d = fo['residual'], fo['row']
        D(f"\\text{{removed equation }} \\phi_{{{d}}}\\text{{: residual }} r_{{{d}}} = \\sum_j K_{{{d}j}}\\,a_j - F_{{{d}}} = {tex(rd) if fo['exact'] else dec(rd)}\\qquad \\text{{(not satisfied automatically)}}")
        yb_tex = tex(fo['yb']) if fo['exact'] and fo['yb'] is not None else dec(fo['yb'])
        scaled = '' if fo['exact'] else f",\\quad \\frac{{|\\mathbf{{y}}^{{\\mathsf{{T}}}}\\mathbf{{b}}|}}{{\\|\\mathbf{{y}}\\|\\,\\|\\mathbf{{b}}\\|}} = {dec(fo['scaled'])}"
        D(f"\\mathbf{{y}}^{{\\mathsf{{T}}}}\\mathbf{{A}} = 0,\\ \\mathbf{{y}}^{{\\mathsf{{T}}}}\\mathbf{{b}} = {yb_tex}{scaled}"
          f"\\qquad (\\mathbf{{A}}\\mathbf{{a}} = \\mathbf{{b}}:\\ {dofs}\\times{dofs - 1}\\ \\text{{full system, compatibility with the left null space}})")
        if fo['compatible']:
            D('\\text{the full system is consistent: the removed equation is also satisfied; the row choice does not change the result}')
        elif fo['compatible'] is False:
            D('\\text{the full system is inconsistent: removing another row (another weighting space) gives another } \\mathbf{a}')
    if form['ivp_slope']:
        sl = form['ivp_slope']
        D(f"u_h'({sp.latex(sl['point'])}) = {tex_value(sl['computed'])},\\qquad u'({sp.latex(sl['point'])}) = {sp.latex(sl['given'])}"
          f"\\qquad \\text{{(the u' condition is satisfied weakly, through the flux)}}")

    # --- yeniden kurulum
    section('Reconstruction of the approximate solution')
    D('u_h^{(e)}(x) = \\sum_{k=1}^{p+1} a^{(e)}_k\\,N_k^{(e)}(x),\\qquad a^{(e)}_k = a_{(e-1)p+k},\\qquad u_h(x) = u_h^{(e)}(x)\\ (x \\in \\Omega_e)')
    for piece in form['reconstruction']:
        e, xe = piece['e'], piece['xe']
        if piece['N'] is not None:
            sub = ' + '.join(f'\\left({tex(sp.sympify(v))}\\right)\\left({sp.latex(Nk)}\\right)' for v, Nk in zip(piece['values'], piece['N']))
            if visible_len(sub + tex(piece['piece'])) <= LINE_LIMIT:
                D(f'u_h^{{({e + 1})}}(x) = {sub} = {tex(piece["piece"])},\\qquad {sp.latex(xe)} \\le x \\le {sp.latex(xe + h)}')
            else:                                         # uzun ifade tek satırda okunmaz kadar küçülmesin: iki satır
                terms = [f'\\left({tex(sp.sympify(v))}\\right)\\left({sp.latex(Nk)}\\right)' for v, Nk in zip(piece['values'], piece['N'])]
                D(f'\\begin{{aligned}} u_h^{{({e + 1})}}(x) &= ' + ' \\\\ &\\quad + '.join(terms)
                  + f' \\\\ &= {tex(piece["piece"])},\\qquad {sp.latex(xe)} \\le x \\le {sp.latex(xe + h)} \\end{{aligned}}')
        else:
            D(f'u_h^{{({e + 1})}}(x) = {tex(piece["piece"])},\\qquad {sp.latex(xe)} \\le x \\le {sp.latex(xe + h)}')
    if n > DETAIL_LIMIT:
        D('\\text{likewise on the other elements}')

    # --- hata ölçüleri
    section('Error measures')
    em = V['errors']
    if em is None:
        D('\\text{not computed: no verified analytical solution (' + ('candidate rejected' if exact['status'] == 'rejected' else exact['reason']) + ')}')
    else:
        D(f'u(x) = {sp.latex(solution)}\\qquad \\text{{({level[exact["status"]]})}}')
        D('\\begin{aligned} &\\max_j |u(x_j) - a_j| &&\\text{nodal error: at the nodes only (may show superconvergence in 1D)} \\\\ '
          '&\\max_x |u - u_h| &&\\text{maximum pointwise error: including element interiors (by sampling; a lower estimate)} \\\\ '
          '&\\|u - u_h\\|_{L^2} &&\\text{average (mean-square) error size} \\\\ '
          '&|u - u_h|_{H^1} &&\\text{error of the derivative (seminorm related to the energy norm)} \\end{aligned}')
        D(',\\quad '.join(f'{name[1]}\\ \\checkmark\\ \\text{{({EVIDENCE_LABEL.get(how, how)})}}' for name, ok, how in exact['checks']))
        exact_or_dec = lambda v: sp.latex(v) if sp.sympify(v).is_Rational else dec(sp.sympify(v).evalf())
        table = ' \\\\ '.join(f'{j + 1} & {tex(xj)} & {exact_or_dec(aj)} & {exact_or_dec(uj)} & {exact_or_dec(ej)}'
                              for j, xj, aj, uj, ej in em['nodal_table'])
        D(f'\\begin{{array}}{{c|c|c|c|c}} j & x_j & a_j = u_h(x_j) & u(x_j) & |u - u_h| \\\\ \\hline {table} \\end{{array}}\\qquad \\text{{(nodal error)}}')
        mids = ' \\\\ '.join(f'{tex(xm)} & {dec(um)} & {dec(uhm)} & {dec(um - uhm)}' for xm, um, uhm in em['midpoints'])
        D(f'\\begin{{array}}{{c|c|c|c}} x & u(x) & u_h(x) & u - u_h \\\\ \\hline {mids} \\end{{array}}\\qquad \\text{{(pointwise error at the midpoint of the first two nodes of each element, }} x_e + h/2p\\text{{)}}')
        rel = lambda v: dec(v) if v is not None else '\\text{undefined}'
        D(f'\\max_j |u(x_j) - a_j| = {dec(em["nodal"])}\\ (\\text{{relative }} {rel(em["rel_nodal"])})')
        D(f'\\max_x |u - u_h| \\approx {dec(em["pointwise"])}\\ (\\text{{relative }} {rel(em["rel_pointwise"])};\\ \\text{{200 samples per element, a lower estimate}})')
        D(f'\\|u - u_h\\|_{{L^2}} = \\left(\\int_{{{X0}}}^{{{X1}}}(u - u_h)^2\\,dx\\right)^{{1/2}} \\approx {dec(em["l2"])}\\ (\\text{{relative }} {rel(em["rel_l2"])})')
        D(f'|u - u_h|_{{H^1}} = \\left(\\int_{{{X0}}}^{{{X1}}}(u\' - u_h\')^2\\,dx\\right)^{{1/2}} \\approx {dec(em["h1"])}\\ (\\text{{relative }} {rel(em["rel_h1"])})')
        if em['nodal'] < 1e-12 * max(1.0, abs(float(max(em['l2'], em['pointwise'])))) and em['pointwise'] > 1e-12:
            D('\\text{no error at the nodes but error inside the elements: accurate nodal values do not mean an accurate function}')
    conv = V['convergence']
    if conv.get('rows'):
        table = ' \\\\ '.join(f"{row['n']} & {dec(row['h'])} & {dec(row['l2'])} & {dec(row['rate_l2']) if row.get('rate_l2') is not None else '-'} & "
                              f"{dec(row['h1'])} & {dec(row['rate_h1']) if row.get('rate_h1') is not None else '-'}" for row in conv['rows'])
        D(f'\\begin{{array}}{{c|c|c|c|c|c}} n & h & \\|e\\|_{{L^2}} & r_{{L^2}} & |e|_{{H^1}} & r_{{H^1}} \\\\ \\hline {table} \\end{{array}}'
          '\\qquad r = \\frac{\\log(E_h/E_{h/2})}{\\log 2}\\ \\text{(independent solver ifem)}')
        expected = conv['expected']
        D('\\text{beklenen: }' + (f"L^2\\ {expected['l2']},\\ H^1\\ {expected['h1']}\\ \\text{{({expected['basis']})}}" if expected
                                 else '\\text{no theoretical rate is claimed for this class}'))
        D(text_block(conv['detail']))
    else:
        D(text_block('convergence: ' + conv['detail']))

    # --- doğrulama ve tanılar
    section('Verification and numerical diagnostics')
    for it in V['items']:
        D(f'{STATUS_MARK[it["status"]]}\\quad \\text{{{it["name"]} ({STATUS_LABEL[it["status"]]})}}')
        D(text_block(it['detail']))                       # uzun ayrıntı tek satırda küçülmesin: satırlara bölünür
    dist = V['distributional']
    if dist and dist.get('jumps') is not None and dist['jumps']:
        shown_j = dist['jumps'][:6]
        D(',\\ '.join(f'J_{{{(k + 1) * degree + 1}}} = {dec(v)}' for k, v in enumerate(shown_j)) + ('\\ \\dots' if len(dist['jumps']) > 6 else '')
          + '\\qquad \\text{(at interior nodes } \\alpha[u_h\'] \\text{; global numbering)}')
        for side, (computed, given) in dist['ends'].items():
            P = X0 if side == 'left' else X1
            D(f"\\alpha u_h'({P}) = {dec(computed)},\\quad Q = {dec(given)}\\qquad \\text{{(the natural condition is satisfied only weakly by the discrete solution)}}")

    # --- uyarılar
    section('Warnings and limitations')
    warnings = []
    g = r['degeneracy']
    if g['kind'] != 'uniform':
        warnings.append('the leading coefficient vanishes in the interval' + (' and changes sign' if g['kind'] == 'sign_changing' else '')
                        + ': no uniform ellipticity/coercivity; convergence theory and stability are not guaranteed')
    for it in V['items']:
        if it['status'] == 'warning':
            warnings.append(f"{it['name']}: {it['detail']}")
    if order == 1:
        warnings.append('first-order standard Galerkin gives no stability guarantee; oscillations may appear in stiff problems')
    if r['ivp']:
        warnings.append("in an initial value problem the u' condition is satisfied weakly; on a coarse mesh an oscillatory solution accumulates phase error")
    if not r['symbolic'] or 'numeric' in methods:
        warnings.append('some integrals were computed numerically (≈); those values have 15 digits')
    unavailable = [it['name'] for it in V['items'] if it['status'] == 'not_run']
    if unavailable:
        warnings.append('checks not run: ' + ', '.join(unavailable))
    if V['errors'] is not None:
        warnings.append('the maximum pointwise error is based on sampling (a lower estimate)')
    for w in warnings or ['no specific warnings']:
        D(text_block(w))
    return '\n\n'.join(blocks)


# ---------- 6. HAT, LATEX BELGESİ VE GİRİŞ ----------
# pdflatex Unicode simgeleri kendiliğinden dizemez; kâğıt metinlerinde (doğrulama ayrıntıları, iz tablosu) geçenler burada
# matematik karşılıklarına bağlanır. MathJax aynı simgeleri doğrudan çizer.
UNICODE_TEX = {'≤': '\\le', '≥': '\\ge', '≠': '\\ne', '≈': '\\approx', '‖': '\\|', 'Σ': '\\Sigma', 'κ': '\\kappa', 'η': '\\eta',
               '∞': '\\infty', 'λ': '\\lambda', 'α': '\\alpha', 'β': '\\beta', 'γ': '\\gamma', 'δ': '\\delta', 'φ': '\\varphi',
               'π': '\\pi', '−': '-', '×': '\\times', '→': '\\to', '·': '\\cdot', '²': '^{2}', '³': '^{3}', '₀': '_{0}', '₁': '_{1}',
               '′': "'", 'ₚ': '_{p}', 'ᵖ': '^{p}', 'ₗ': '_{l}', 'ₑ': '_{e}', 'ₖ': '_{k}', 'ᵢ': '_{i}', 'ⱼ': '_{j}',
               **{chr(0x2082 + k): f'_{{{k + 2}}}' for k in range(8)}, '≡': '\\equiv', '≢': '\\not\\equiv', '⇒': '\\Rightarrow', '∫': '\\int', 'χ': '\\chi', 'Ω': '\\Omega'}


def build_tex(text):
    body = []
    for block in text.split('\n\n'):
        heading = re.match(r'#+\s+(.*)', block)
        display = re.fullmatch(r'\$\$(.*)\$\$', block, re.S)
        body.append(f'\\section*{{{heading.group(1)}}}' if heading else f'\\[ \\fit{{{display.group(1)}}} \\]' if display else block)
    document = '\n\n'.join(body)
    # Kullanıcı girdisi \texttt içinde aynen durur; ^ metin kipinde LaTeX hatasıdır, karakter olarak basılır.
    document = re.sub(r'\\texttt\{([^{}]*)\}', lambda m: '\\texttt{' + m.group(1).replace('^', '\\textasciicircum{}') + '}', document)
    return '\n'.join([
        '\\documentclass[11pt,a4paper]{article}',
        '\\usepackage[utf8]{inputenc}', '\\usepackage[T1]{fontenc}', '\\usepackage{lmodern,amsmath,amssymb,graphicx}',
        '\\usepackage[margin=18mm]{geometry}', '\\setlength{\\parindent}{0pt}', '\\setcounter{MaxMatrixCols}{20}',
        *[f'\\DeclareUnicodeCharacter{{{ord(ch):04X}}}{{\\ensuremath{{{tex_}}}}}' for ch, tex_ in UNICODE_TEX.items()],
        '\\newcommand{\\fit}[1]{\\resizebox{\\ifdim\\width>\\linewidth\\linewidth\\else\\width\\fi}{!}{$\\displaystyle #1$}}',
        '\\begin{document}', '{\\LARGE\\textbf{FEM LAB --- Finite Element Solution Sheet}}\\par\\smallskip\\hrule',
        document, '\\end{document}'])


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


# Katman 2'nin tarayıcıya giden kısmı: iz tablosu ve katsayıların 9 noktadaki değerleri. Tarayıcı kendi (üçüncü, bağımsız)
# ayrıştırıcısıyla aynı noktaları hesaplayıp karşılaştırır: üç ayrıştırıcı aynı problemi okumuş olmalı.
def trace_payload(valid):
    P = valid['problem']
    points = [P['x0'] + (P['x1'] - P['x0']) * sp.Rational(2 * k + 1, 18) for k in range(9)]
    samples = {key: [float(P[key].subs(x, pt).evalf(20)) for pt in points] for key, _ in COEFFICIENTS}
    return {'trace': [{k: v for k, v in row.items()} for row in valid['trace']],
            'coefficient_samples': {'points': [float(pt) for pt in points], **samples}}


def verification_payload(V):
    keep = lambda it: {k: v for k, v in it.items() if k in ('name', 'status', 'detail')}
    return [keep(it) for it in V['items']]


# Katman 3'ü çalıştırır; motorun ayrık teşhisi (ör. tekil sistem) gelirse Katman 4'ün sürekli sınıflandırması yine çalışır:
# sürekli problem sembolik olarak çözümsüz/teksiz kanıtlanmışsa kullanıcıya o (daha temel) neden bildirilir; değilse ayrık
# teşhis, sürekli bilgiyle birlikte yükselir. Ayrık tekillik tek başına "diferansiyel problem çözümsüz" diye sunulmaz.
def engine_or_diagnosis(valid, symbolic):
    import fem_verify
    try:
        return fem_engine(valid, symbolic)
    except Diagnosis as error:
        P = valid['problem']
        c = fem_verify.classify_problem(P['alpha'], P['beta'], P['gamma'], P['f'], P['x0'], P['x1'], P['conditions'],
                                        P['slopes'], P['order'], P['ivp'], P['degeneracy'])
        stop = fem_verify.continuous_diagnosis(c)
        if stop is not None:
            raise stop
        error.classification, error.degeneracy = c, P['degeneracy']
        raise


# Hat (pipeline): Katman 1 → 2 → 3 → 4. Bir katman başarısızsa sonraki çalışmaz (UserError / Diagnosis yükselir).
def run_pipeline(data, symbolic=True):
    import fem_verify
    valid = validate_problem(define_problem(data))
    r = engine_or_diagnosis(valid, symbolic)
    V = fem_verify.verify(r)
    return valid, r, V


def make_sheet(data, symbolic=True):
    valid, r, V = run_pipeline(data, symbolic)
    text = build_sheet(valid, r, V)
    solution = V['exact']['u']
    # K: tarayıcının montajladığı matrisle girdi girdi karşılaştırılsın diye (büyük sistemde yük olmasın: ≤ 65 bilinmeyen).
    K = [[float(v) for v in row] for row in r['K'].evalf(17).tolist()] if r['dofs'] <= 65 else None
    # Eleman matrisleri ve yük vektörleri (motorun gerçekten montajladığı): bağımsız eleman düzeyi karşılaştırma için.
    element_K = [[[float(v) for v in row] for row in rec['K'].evalf(17).tolist()] for rec in r['elements']]
    element_F = [[float(v) for v in rec['F'].evalf(17)] for rec in r['elements']]
    return {'text': text, 'tex': build_tex(text), 'nodal': [float(sp.sympify(v).evalf()) for v in r['a']],
            'exact': str(solution) if solution is not None else None,
            'exact_js': exact_for_browser(solution, r['x0'], r['x1']), 'K': K, 'element_K': element_K, 'element_F': element_F,
            **verdict_payload(V['exact']),
            'fixed': {str(i): float(v) for i, v in r['fixed'].items()},
            'first_order_residual': float(r['first_order']['residual']) if r['first_order'] else None,
            **classification_payload(V['classification'], r['degeneracy']), **trace_payload(valid),
            'verification': verification_payload(V)}


# Katman 1–3'ü (classify=True ise sürekli sınıflandırmayı da) çalıştıran kısa yol: testler ve hızlı analitik istek için.
# classify=True iken sembolik kanıtlı sürekli K1/K2, hat ile aynı biçimde Diagnosis olarak yükselir.
def solve_problem(data, symbolic=True, classify=False):
    valid = validate_problem(define_problem(data))
    if not classify:
        return fem_engine(valid, symbolic)
    import fem_verify
    r = engine_or_diagnosis(valid, symbolic)
    r['classification'] = fem_verify.classify(r)
    stop = fem_verify.continuous_diagnosis(r['classification'])
    if stop is not None:
        raise stop
    return r


def exact_only(data):
    """Analitik çözüm + sürekli sınıflandırma (kâğıt yok): grafikler ve teşhis kutusu için hızlı istek. Ağ n = 1, p = 1
    ile en ucuz haliyle kurulur. Sürekli düzeyde kanıtlanmış ret (K1/K2) 'continuous_error' olarak döner; n = 1 ağının
    ayrık tekilliği sürekli problem hakkında bilgi vermediği için sınıflandırmayı silmez."""
    import fem_verify
    try:
        r = solve_problem({**data, 'n': 1, 'degree': 1}, symbolic=False, classify=True)
    except Diagnosis as error:
        if error.level == 'continuous':
            return {'exact_js': None, 'continuous_error': {'category': error.category, 'message': str(error)}}
        return {'exact_js': None, **classification_payload(getattr(error, 'classification', None), getattr(error, 'degeneracy', None))}
    except UserError:
        return {'exact_js': None}
    verdict = fem_verify.exact_solution(r)
    return {'exact_js': exact_for_browser(verdict['u'], r['x0'], r['x1']), **verdict_payload(verdict),
            **classification_payload(r['classification'], r['degeneracy'])}


def main():
    payload = json.load(sys.stdin)
    try:
        if payload.get('mode') == 'exact':
            output, code = exact_only(payload), 0
        else:
            output, code = make_sheet(payload, payload.get('symbolic', True)), 0
    except UserError as error:
        # Teşhisse makine alanları da gider (kategori, düzey): istemci ve korpus çalıştırıcısı mesaj metnini ayrıştırmaz.
        output = {'error': str(error), **({'category': error.category, 'level': error.level} if isinstance(error, Diagnosis) else {})}
        code = 1
    sys.stdout.write(json.dumps(output, ensure_ascii=False))
    sys.stdout.flush()
    os._exit(code)      # süre sınırını aşmış bir dsolve iş parçacığı varsa beklemeden çık


if __name__ == '__main__':
    import fem_sheet    # tek modül kopyası: fem_verify 'fem_sheet'i içe aktarır; __main__ ile ayrı kopya olmasın (Diagnosis aynı sınıf)
    fem_sheet.main()
