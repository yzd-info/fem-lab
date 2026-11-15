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
            return True, 'sembolik'
    except Exception:
        pass
    with mpmath.workdps(30):
        for k in range(9):
            pt = x0 + (x1 - x0) * sp.Rational(2 * k + 1, 18)
            try:
                a, b = complex(sp.N(first.subs(x, pt), 30)), complex(sp.N(second.subs(x, pt), 30))
            except Exception:
                return False, 'değerlendirilemedi'
            if not abs(a - b) <= EQUIVALENCE_TOL * max(1.0, abs(a), abs(b)):
                return False, f'x = {dec(pt)}: {dec(a.real)} ≠ {dec(b.real)}'
    return True, 'sayısal (9 nokta, 30 basamak)'


# Önde gelen katsayının (2. mertebede α, 1. mertebede β) sıfırları. Neden: α'nın sıfırı ya da işaret değişimi standart
# H¹ kuramını (düzgün eliptiklik |α| ≥ α₀ > 0) bozar; ama TEK BAŞINA problemin çözümsüz olduğunu göstermez. Bu yüzden sonuç
# yalnız sınıflandırma ve uyarıdır, ret değildir. Sıfırlar mümkünse sembolik (solveset), değilse sayısal bulunur ve
# hangi yolun kullanıldığı kaydedilir. Sayısal yol çok dar bir değme sıfırını kaçırabilir (bilinen sınır).
# Katman 2'de durur: problemin bir ÖZELLİĞİDİR (sınıfını belirler), çözümün varlığı hakkında hüküm vermez.
def coefficient_degeneracy(lead, x0, x1):
    span = x1 - x0
    if not lead.free_symbols:
        return {'kind': 'uniform', 'points': [], 'endpoints': [], 'interior': [], 'sign_change': [], 'method': 'sembolik'}
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


# Katman 2. Amaç: çözücüye giden problemin kullanıcının yazdığıyla AYNI olduğunu göstermek. Çözümün varlığı/tekliği burada
# SORULMAZ (o Katman 4'ün işi). Her metin alanı iki bağımsız ayrıştırıcıdan geçer (parse: simge simge çevirici;
# independent_parse: SymPy dilbilgisi) ve sonuç iz tablosuna yazılır. Kural ihlali ya da eşdeğerliğin kurulamaması hesabı
# DURDURUR (UserError). Girdi asla sessizce düzeltilmez; büyük/küçük harf birleştirmesi tek dönüşümdür ve izde yazılır.
def validate_problem(raw):
    trace = []

    def row(item, given, internal, ok, how):
        trace.append({'item': item, 'raw': given, 'internal': internal, 'ok': bool(ok), 'how': how})
        if not ok:
            raise UserError(f'Problem doğrulaması: {item} girdisi ile iç gösterim eşdeğer değil ({how}). Hesap durduruldu.')

    def second(text, label):
        try:
            return independent_parse(text)
        except Exception:
            raise UserError(f'Problem doğrulaması: {label} ikinci ayrıştırıcıyla okunamadı; eşdeğerlik kurulamadı. Hesap durduruldu.')

    shown = lambda text: text if text.strip() else ''
    x0, x1 = parse_number(raw['x0'], 'x₀'), parse_number(raw['x1'], 'x₁')
    for label, text, value in (('x₀', raw['x0'], x0), ('x₁', raw['x1'], x1)):
        ok, how = same_expression(value, second(text, label), 0, 1)
        row(label, shown(text), sp.latex(value), ok, how)
    if not x0 < x1:
        raise UserError('Aralıkta x₀, x₁\'den küçük olmalı.')

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
            raise UserError('n ve p tam sayı olmalı.')
        if not number.is_integer():
            raise UserError('n ve p tam sayı olmalı.')
        return int(number)
    n, degree = whole(raw['n']), whole(raw['degree'])
    if not (1 <= n <= 64 and 1 <= degree <= 4):
        raise UserError('n 1–64, p 1–4 aralığında olmalı.')
    order = 2 if alpha != 0 else 1 if beta != 0 else 0
    if order == 0:
        raise UserError('α ve β birlikte sıfır olamaz: bu bir diferansiyel denklem değil.')

    # Koşul türleri: none · dirichlet (u = g) · neumann (u' = g, TÜREV koşulu) · flux (α u' = q, AKI koşulu) ·
    # both (aynı uçta u ve u': başlangıç/son değer problemi; u' değeri slope alanında). Türev ve akı koşulu birbirine
    # dönüştürülmez: zayıf forma giren akı neumann'da Q = α(x_b)·g, flux'ta Q = q'dur ve iz tablosunda açıkça yazılır.
    conditions, slopes = {}, {}
    for side, point in (('left', x0), ('right', x1)):
        given = raw[side]
        kind = given['type']
        if kind not in CONDITION_KINDS:
            raise UserError('Bilinmeyen koşul türü.')
        value = parse_number(given['value'], 'koşul değeri') if kind != 'none' else None
        conditions[side] = (kind, value)
        if kind == 'both':
            slopes[side] = parse_number(given['slope'], "u' koşul değeri")
    kinds = [conditions['left'][0], conditions['right'][0]]
    if order == 1 and (any(k in kinds for k in ('neumann', 'flux', 'both')) or kinds.count('dirichlet') != 1):
        raise UserError('Birinci mertebe denklem tam olarak bir koşul alır: uçlardan birinde u değeri ver, diğerini boş bırak.')
    # Başlangıç değer problemi: iki koşul aynı uçta. Öbür uca koşul verilirse problem aşırı belirlenir (3 koşul, 2. mertebe).
    ivp = 'both' in kinds
    if ivp and sorted(kinds) != ['both', 'none']:
        raise UserError("Bir uçta hem u hem u' verildiyse (başlangıç değer problemi) diğer uca koşul verilmez.")
    # u' (ya da akı) koşulu zayıf forma yalnızca akı Q = α u' üzerinden girer. α(uç) = 0 ise verilen değer sessizce
    # kaybolurdu (u' koşulunda) ya da 0 = q çelişkisi olurdu (akı koşulunda). Bu yüzden o uçta α ≠ 0 şartı aranır.
    for side, point in (('left', x0), ('right', x1)):
        if conditions[side][0] in ('neumann', 'both') and alpha.subs(x, point) == 0:
            raise UserError("α bu uçta sıfır: u' koşulu denkleme girmez (α u' akısı sıfır olur).")
        if conditions[side][0] == 'flux' and alpha.subs(x, point) == 0:
            raise UserError("α bu uçta sıfır: akı koşulu α u' = q ancak q = 0 iken anlamlıdır ve o zaman doğal koşulun kendisidir; "
                            "bu uçta koşulu 'verilmedi' bırak ya da α'yı kontrol et.")

    for side, point, name in (('left', x0, 'x₀ koşulu'), ('right', x1, 'x₁ koşulu')):
        kind, value = conditions[side]
        given, P = raw[side], sp.latex(point)
        a_b = alpha.subs(x, point)
        if kind == 'none':
            meaning = ('\\text{koşul yok (1. mertebe)}' if order == 1 else '\\text{serbest uç (başlangıç değer problemi)}' if ivp
                       else '\\alpha(' + P + ') = 0\\ \\text{(dejenere uç; akı kendiliğinden sıfır)}' if a_b == 0
                       else f"\\alpha u'({P}) = 0\\ \\text{{(doğal koşul)}}")
            trace.append({'item': name, 'raw': '(verilmedi)', 'internal': meaning, 'ok': True, 'how': 'tür: verilmedi'})
            continue
        ok, how = same_expression(value, second(given['value'], name), x0, x1)
        if kind == 'dirichlet':
            internal = f'u({P}) = {sp.latex(value)}\\ \\text{{(temel; akı bilinmiyor)}}'
        elif kind == 'neumann':
            internal = f"u'({P}) = {sp.latex(value)}\\ \\Rightarrow\\ Q = \\alpha({P})\\,u'({P}) = {sp.latex(sp.simplify(a_b * value))}"
        elif kind == 'flux':
            internal = f"\\alpha({P})\\,u'({P}) = {sp.latex(value)}\\ \\Rightarrow\\ Q = {sp.latex(value)}"
        else:
            ok2, how2 = same_expression(slopes[side], second(given['slope'], name), x0, x1)
            ok, how = ok and ok2, how if ok2 else how2
            internal = (f"u({P}) = {sp.latex(value)},\\ u'({P}) = {sp.latex(slopes[side])}\\ \\Rightarrow\\ "
                        f"Q = {sp.latex(sp.simplify(a_b * slopes[side]))}")
        label = {'dirichlet': 'u =', 'neumann': "u' =", 'flux': "α u' =", 'both': "u ve u' ="}[kind]
        row(name, f"{label} {given['value']}" + (f", {given['slope']}" if kind == 'both' else ''), internal, ok, how)

    present = ', '.join(f'{name} {"≢" if coef[key] != 0 else "≡"} 0' for key, name in (('alpha', 'α'), ('beta', 'β'), ('gamma', 'γ'), ('f', 'f')))
    trace.append({'item': 'operatör', 'raw': f"(({shown(raw['alpha']) or '0'}) u')' + ({shown(raw['beta']) or '0'}) u' + "
                  f"({shown(raw['gamma']) or '0'}) u = {shown(raw['f']) or '0'}",
                  'internal': tex_equation(alpha, beta, gamma, f), 'ok': True, 'how': f'{present}; mertebe {order}'})
    trace.append({'item': 'ağ', 'raw': f"n = {raw['n']}, p = {raw['degree']}", 'internal': f'n = {n},\\ p = {degree}',
                  'ok': True, 'how': 'tam sayı'})

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
    return sp.Symbol(f'\\hat{{a}}^{{(e)}}_{{{k}}}')      # e elemanının yerel katsayısı (global a_j ile karıştırılmasın)


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
            raise RuntimeError(f'Türetilen N_{k + 1} hesapta kullanılanla aynı değil')
    delta = sp.Matrix(degree + 1, degree + 1, lambda i, j: shapes[i].subs(s, nodes[j]))
    partition = sp.simplify(sum(shapes))
    completeness = [sp.simplify(sum(Nk * nodes[k] ** m for k, Nk in enumerate(shapes)) - s ** m) == 0 for m in range(degree + 1)]
    if delta != sp.eye(degree + 1) or partition != 1 or not all(completeness):
        raise RuntimeError('Şekil fonksiyonları δ_ij / birim bölüşüm / tamlık özelliğini sağlamıyor')
    return dict(degree=degree, h=h, cs=cs, hats=hats, nodes=nodes, poly=poly, conds=conds, V=V, sol=sol,
                W=V[1:, 1:], rhs=[hv - hats[0] for hv in hats[1:]], shapes=shapes,
                derivatives=[sp.diff(Nk, s) for Nk in shapes], delta=delta, partition=partition, completeness=completeness)


# Tek bir eleman integrali. Dönüş YAPILANDIRILMIŞ kayıttır: çarpanlar, açılmış integrand, ilkel fonksiyon, değer, yöntem
# ('sıfır' / 'polinom' / 'sembolik' / 'sayısal'). Kapalı biçim bulunamazsa mpmath.quad (uyarlamalı tanh–sinh) ile 15 basamak.
def integral(label, terms, h, symbolic):
    integrand = sum(sign * sp.prod(factors) for sign, factors in terms)
    record = {'label': label, 'terms': terms, 'h': h, 'integrand': None, 'primitive': None}
    polynomial = integrand.is_polynomial(s)
    if polynomial:
        integrand = sp.expand(integrand)
    if integrand == 0:
        return {**record, 'method': 'sıfır', 'value': sp.Integer(0)}
    if symbolic:
        try:
            if polynomial:
                primitive = sp.integrate(integrand, s)
                value = sp.simplify(primitive.subs(s, h) - primitive.subs(s, 0))
                return {**record, 'method': 'polinom', 'integrand': integrand, 'primitive': primitive, 'value': value}
            value = sp.simplify(sp.integrate(integrand, (s, 0, h)))
            if not value.has(sp.Integral) and value.is_finite and value.is_real:
                return {**record, 'method': 'sembolik', 'value': value}
        except Exception:
            pass
    value = sp.Float(mpmath.quad(sp.lambdify(s, integrand, 'mpmath'), [0, float(h)]), 15)
    return {**record, 'method': 'sayısal', 'value': value}


# Tam (rasyonel/sembolik) yolda tekil sistemin sınıflandırılması: rank ve uyumluluk KESİN (tolerans yok).
def classify_discrete_exact(K_ff, F_f, ivp, has_dirichlet):
    if ivp:
        return Diagnosis('K3', 'ayrık sistem', IVP_COARSE_DETAIL)
    rank, rank_aug, size = K_ff.rank(simplify=True), K_ff.row_join(F_f).rank(simplify=True), K_ff.shape[1]
    if rank_aug > rank:
        example = '' if has_dirichlet else " Örn. saf Neumann'da ∫f ≠ akı farkı."     # örnek yalnız u hiç verilmemişse ilgili
        return Diagnosis('K1', 'ayrık sistem', f'rank K = {rank} < rank [K|F] = {rank_aug} (tam aritmetik): sağ taraf uyumsuz. '
                                                f'Bu AYRIK sistemin hükmüdür; sürekli problem hakkında ayrıca bkz. sınıflandırma.{example}')
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
    exact = not (K_ff.has(sp.Float) or F_f.has(sp.Float)) and len(free) <= EXACT_SOLVE_LIMIT
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
             exact=exact, reactions=reactions, reaction_slopes=reaction_slopes, symbolic=symbolic, condition=condition, degeneracy=degeneracy, classification=None)
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
                raise RuntimeError(f'Galerkin denklemi i={i + 1}, K satırıyla tutmuyor')
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
STATUS_MARK = {'geçti': '\\checkmark', 'başarısız': '\\times', 'uyarı': '!', 'yapılmadı': '-'}


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
    if step['method'] == 'sıfır':
        return wrap('= 0')
    if step['method'] == 'polinom':
        return wrap(f"= \\int_{{0}}^{{{upper}}}\\left({sp.latex(step['integrand'])}\\right)ds = "
                    f"\\left[{sp.latex(step['primitive'])}\\right]_{{0}}^{{{upper}}} = {tex_value(step['value'])}")
    if step['method'] == 'sembolik':
        return wrap(f"= {tex_value(step['value'])}")
    return wrap(f"\\approx {dec(step['value'])}")


# Motorun şekil fonksiyonu kurulum kaydını (Vandermonde yolu) blok blok yazar.
def format_shape(c):
    degree, h, cs, hats = c['degree'], c['h'], c['cs'], c['hats']
    B = []
    B.append(f'u^{{e}}(s) = {sp.latex(c["poly"])},\\qquad 0 \\le s \\le h = {sp.latex(h)}\\qquad '
             f'\\text{{({degree}. derece polinom: {degree + 1} katsayı, {degree + 1} düğüm)}}')
    B.append('s_k = (k - 1)\\,\\frac{h}{p}:\\qquad ' + ',\\ '.join(f's_{{{k + 1}}} = {sp.latex(v)}' for k, v in enumerate(c['nodes'])))
    B.append('\\begin{aligned} ' + ' \\\\ '.join(f'u^{{e}}(s_{{{k + 1}}}) &= {sp.latex(cond)} = {sp.latex(hats[k])}'
                                               for k, cond in enumerate(c['conds'])) + ' \\end{aligned}\\qquad \\text{(interpolasyon koşulları)}')
    B.append(f'{tex_matrix(c["V"].tolist())}\\,{tex_column(cs)} = {tex_column(hats)}\\qquad \\mathbf{{V}}\\mathbf{{c}} = \\hat{{\\mathbf{{a}}}}^{{(e)}}'
             f'\\quad \\text{{(interpolasyon matrisi }} \\mathbf{{V}}\\text{{, }} V_{{km}} = s_k^{{\\,m}}\\text{{)}}')
    B.append(f's_1 = 0 \\;\\Rightarrow\\; c_0 = {sp.latex(c["sol"][cs[0]])}')
    B.append(f'{tex_matrix(c["W"].tolist())}\\,{tex_column(cs[1:])} = {tex_column(c["rhs"])}\\qquad \\text{{(kalan {degree} denklem, }} c_0 \\text{{ yerine konmuş)}}')
    B.append('\\begin{aligned} ' + ' \\\\ '.join(f'{sp.latex(cm)} &= {sp.latex(sp.collect(sp.expand(c["sol"][cm]), hats))}' for cm in cs)
             + ' \\end{aligned}\\qquad \\text{(çözüm)}')
    B.append('u^{e}(s) = ' + ' + '.join(f'{sp.latex(hv)}\\left({sp.latex(Nk)}\\right)' for hv, Nk in zip(hats, c['shapes']))
             + f' = \\sum_{{k=1}}^{{{degree + 1}}} \\hat{{a}}^{{(e)}}_k\\,N_k(s)\\qquad \\text{{(}}\\hat{{a}}^{{(e)}}_k\\text{{ katsayılarına göre gruplandı)}}')
    B.append('\\begin{aligned} ' + ' \\\\ '.join(f"N_{{{k + 1}}}(s) &= {sp.latex(Nk)}, &\\qquad N_{{{k + 1}}}'(s) &= {sp.latex(dNk)}"
                                               for k, (Nk, dNk) in enumerate(zip(c['shapes'], c['derivatives']))) + ' \\end{aligned}')
    B.append(f'\\big[N_i(s_j)\\big] = {tex_matrix(c["delta"].tolist())} = \\mathbf{{I}}\\quad (N_i(s_j) = \\delta_{{ij}}),\\qquad '
             f'\\sum_k N_k(s) = {sp.latex(c["partition"])}\\ \\text{{(birim bölüşüm)}}')
    B.append(f'\\sum_k N_k(s)\\,s_k^{{\\,m}} = s^m,\\quad m = 0, \\dots, {degree}\\qquad \\text{{(polinom tamlığı: {degree}. dereceye kadar her polinom tam temsil edilir)}}')
    return B


# Testlerin ve dışarıdan çağıranların kolaylığı için: kurulumu motordan alıp biçimler (kâğıt kendisi r['shape']'i kullanır).
def shape_function_derivation(degree, h, N_engine):
    return format_shape(shape_function_construction(degree, h, N_engine))


# Sürekli sınıflandırma (Katman 4) ve katsayı dejenereliği (Katman 2) kanıt düzeyiyle yazılır.
def classification_lines(r, c):
    lines, g = [], r.get('degeneracy')
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


def raw_tex(text):
    """Kullanıcının yazdığı metin, karakteri karakterine (\\texttt). LaTeX'e özel karakterler kaçırılır; ^ build_tex'te çevrilir."""
    if not text:
        return '\\text{(boş)}'
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
            return f'u({P}) = {sp.latex(value)}\\ \\text{{(temel koşul)}}'
        if kind == 'neumann':
            return f"u'({P}) = {sp.latex(value)}\\ \\text{{(türev koşulu; akı }} Q = \\alpha({P})\\,u'({P})\\text{{)}}"
        if kind == 'flux':
            return f"\\alpha({P})\\,u'({P}) = {sp.latex(value)}\\ \\text{{(akı koşulu; }} Q = {sp.latex(value)}\\text{{)}}"
        if kind == 'both':
            return f"u({P}) = {sp.latex(value)},\\ u'({P}) = {sp.latex(r['slopes'][side])}\\ \\text{{(başlangıç koşulları)}}"
        if r['ivp'] or order != 2:
            return None
        if alpha.subs(x, point) == 0:
            return f"\\alpha({P}) = 0\\ \\text{{(verilmedi: akı terimi kendiliğinden sıfır, bağımsız koşul değil)}}"
        return f"\\alpha u'\\big|_{{x={P}}} = 0\\ \\text{{(verilmedi: doğal koşul)}}"

    exact = V['exact']
    solution = exact['u']
    level = {'verified_symbolic': 'sembolik doğrulandı', 'verified_numeric': 'sayısal doğrulandı (50 basamak; kanıt değil)'}

    # --- Verilenler
    lhs_tex = tex_equation(alpha, beta, gamma, f).split(' = ')[0]
    add('## Verilenler')
    mark = {True: '\\checkmark', False: '\\times'}
    rows = ' \\\\ '.join(f"\\text{{{t['item']}}} & {raw_tex(t['raw'])} & {t['internal']} & {mark[t['ok']]}\\ \\text{{{t['how']}}}"
                         for t in valid['trace'])
    D(f'\\begin{{array}}{{l|l|l|l}} \\text{{alan}} & \\text{{girdi (aynen)}} & \\text{{iç gösterim}} & \\text{{eşdeğerlik}} \\\\ \\hline {rows} \\end{{array}}')
    D(tex_equation(alpha, beta, gamma, f) + f', \\qquad x \\in \\Omega = ({X0},\\, {X1})')
    D(f'\\mathcal{{L}}(u) = {lhs_tex}\\qquad \\text{{({order}. mertebe diferansiyel operatör; güçlü form: }} \\mathcal{{L}}(u) = f\\text{{)}}')
    D(',\\qquad '.join(c for c in (condition_tex('left', x0), condition_tex('right', x1)) if c))
    D(f'n = {n},\\qquad h = \\frac{{x_1 - x_0}}{{n}} = {sp.latex(h)},\\qquad p = {degree}\\ \\text{{({DEGREE_NAMES[degree]} eleman: eleman başına {degree + 1} düğüm)}},'
      f'\\qquad N = np + 1 = {dofs}\\ \\text{{(toplam düğüm = bilinmeyen sayısı)}}')
    method_text = ('Bubnov–Galerkin' if not r['ivp'] else 'Petrov–Galerkin (test fonksiyonları akısı bilinmeyen uçta sıfır)')
    D(f'\\text{{Yöntem: {method_text}}}')
    for line in classification_lines(r, V['classification']):
        D(line)
    if solution is not None:
        D(f'u(x) = {sp.latex(solution)}\\qquad \\text{{(analitik çözüm; {level[exact["status"]]}: denklem ve tüm koşullar)}}')
    elif exact['status'] == 'rejected':
        failed = next(c for c in exact['checks'] if not c[1])
        D(f"\\text{{analitik aday }} {sp.latex(exact['candidate'])}\\ \\text{{reddedildi: }} {failed[0][1]}\\ \\text{{sağlanmıyor (kullanılmadı)}}")
    else:
        D(f'\\text{{analitik çözüm yok: {exact["reason"]}}}')

    section('Gösterim')
    D('\\begin{aligned} u(x)&:\\ \\text{tam çözüm (biliniyorsa)} & u_h(x)&:\\ \\text{sonlu eleman yaklaşımı} \\\\ '
      '\\phi_j(x)&:\\ \\text{global yaklaşım (şekil) fonksiyonu} & N_k(s),\\ N_k^{(e)}(x)&:\\ \\text{yerel şekil fonksiyonu, } e \\text{ elemanında} \\\\ '
      'a_j&:\\ \\text{bilinmeyen global düğüm katsayısı} & \\hat{a}^{(e)}_k&:\\ e \\text{ elemanının yerel katsayısı} \\\\ '
      'c_m&:\\ \\text{eleman polinomunun katsayısı} & R(u_h, x)&:\\ \\text{kalıntı} \\\\ '
      'w_i(x)&:\\ \\text{ağırlık (test) fonksiyonu} & K_{ij},\\ F_i&:\\ \\text{sistem matrisi ve yük vektörü girdileri} \\end{aligned}')
    D('\\phi_j(x) = \\begin{cases} N_k^{(e)}(x), & x \\in \\Omega_e,\\ j = (e-1)p + k \\\\ 0, & \\text{diğer elemanlarda} \\end{cases}'
      '\\qquad \\text{(genel Galerkin tabanı } \\phi_j \\text{, sonlu eleman interpolasyonunda parça parça } N_k^{(e)} \\text{ olur)}')

    section('Ağ')
    nodes = form['nodes']
    D('\\big[x_j\\big] = ' + tex_row(nodes) + f'\\qquad \\Omega_e = [x_e,\\ x_e + h],\\ e = 1, \\dots, {n}')
    shown_e = list(range(n)) if n <= 6 else [0, 1, None, n - 1]
    table = ' \\\\ '.join('\\vdots & \\vdots & \\vdots' if e is None else
                          f'{e + 1} & [{sp.latex(x0 + e * h)},\\ {sp.latex(x0 + (e + 1) * h)}] & '
                          + ',\\ '.join(f'{e * degree + k + 1}' for k in range(degree + 1)) for e in shown_e)
    D(f'\\begin{{array}}{{c|c|c}} e & \\Omega_e & \\text{{yerel }} k = 1, \\dots, {degree + 1} \\;\\mapsto\\; \\text{{global }} j = (e-1)p + k \\\\ \\hline {table} \\end{{array}}')

    section('Şekil fonksiyonları (polinom interpolasyonundan)')
    D('s = x - x_e \\in [0, h]\\qquad \\text{(yerel koordinat: her elemanda aynı türetme)}')
    for block in format_shape(r['shape']):
        D(block)
    D('x = x_e + s,\\quad dx = ds,\\quad \\frac{dN_k^{(e)}}{dx} = \\frac{dN_k}{ds}\\ (\\text{Jacobian } 1),\\qquad N_k^{(e)}(x) = N_k(x - x_e),'
      '\\qquad \\hat{a}^{(e)}_k = a_{(e-1)p+k}\\qquad \\text{(yerel → global)}')
    D('\\text{süreklilik: komşu elemanlar ortak düğümde aynı } a_j \\text{ katsayısını kullanır} \\Rightarrow u_h \\in C^0;\\ '
      'u_h\' \\text{ düğümlerde genel olarak süreksiz}')
    for e in form['shown']:
        E = form['elements'][e]
        D(f'e = {e + 1}:\\qquad ' + ',\\quad '.join(f'N_{{{k + 1}}}^{{({e + 1})}}(x) = {sp.latex(Nk)}' for k, Nk in enumerate(E['N']))
          + ',\\qquad ' + ',\\ '.join(f'\\hat{{a}}^{{({e + 1})}}_{{{k + 1}}} = {sp.latex(a)}' for k, a in enumerate(E['coefficients'])))

    section('Yaklaşık çözüm')
    D(f'u_h(x) = \\sum_{{j=1}}^{{{dofs}}} a_j\\,\\phi_j(x)\\qquad \\text{{(}} a_j \\text{{: belirlenecek bilinmeyenler)}}')
    D('\\phi_j(x_i) = \\delta_{ij} \\;\\Rightarrow\\; u_h(x_i) = a_i\\qquad \\text{(katsayılar düğüm değerleridir)}')
    D('x \\in \\Omega_e:\\quad u_h(x) = u_h^{(e)}(x) = \\sum_{k=1}^{p+1} \\hat{a}^{(e)}_k\\,N_k^{(e)}(x)')
    for e in form['shown']:
        E = form['elements'][e]
        D(f'u_h^{{({e + 1})}}(x) = ' + ' + '.join(f'{sp.latex(a)}\\left({sp.latex(Nk)}\\right)' for a, Nk in zip(E['coefficients'], E['N']))
          + f',\\qquad x \\in [{sp.latex(E["xe"])}, {sp.latex(E["xe"] + h)}]')
    if n > len(form['shown']):
        D('\\text{diğer elemanlar aynı yapıda}')

    # --- Galerkin ağırlıklı kalıntı
    section('Galerkin Weighted-Residual Formulation (Galerkin ağırlıklı kalıntı)')
    D(f'R(u_h, x) = \\mathcal{{L}}(u_h) - f = \\left[{lhs_tex}\\right]_{{u = u_h}}' + ('' if f == 0 else f' - \\left({sp.latex(f)}\\right)'))
    interior_name = 'R_e(x)' if order == 2 else 'R(u_h, x)\\big|_{\\Omega_{e}}'
    for e in form['shown']:
        label = f'R_{{{e + 1}}}(x)' if order == 2 else f'R(u_h, x)\\big|_{{\\Omega_{{{e + 1}}}}}'
        D(f'{label} = {sp.latex(form["elements"][e]["residual"])}' + ('\\qquad \\text{(eleman içi klasik kalıntı)}' if order == 2 else ''))
    if order == 2:
        D("\\mathcal{L}(u_h) - f = \\sum_e R_e\\,\\chi_{\\Omega_e} + \\sum_{k} J_k\\,\\delta(x - x_k),\\qquad "
          "J_k = \\alpha(x_k)\\big[u_h'(x_k^+) - u_h'(x_k^-)\\big]\\qquad \\text{(dağılımsal kalıntı; iç düğümlerde türev sıçraması)}")
        D('\\text{eleman içi } R_e \\text{ tek başına global kalıntı değildir; } C^0 \\text{ elemanda düğüm terimleri } J_k \\text{ de vardır}')
    D(f'\\int_{{{X0}}}^{{{X1}}} w_i(x)\\,R(u_h, x)\\,dx = 0\\qquad \\text{{(kalıntı, ağırlık fonksiyonlarına göre ortalamada sıfır)}}')
    D('w_i(x) = \\phi_i(x)\\qquad \\text{(' + ('Bubnov–Galerkin: ağırlık fonksiyonları yaklaşım fonksiyonlarıyla aynı' if not r['ivp']
      else 'Petrov–Galerkin: ağırlıklar yaklaşım fonksiyonlarından seçilir, test uzayı farklı') + ')}')
    tests_used = list(r['rows'])
    D(f'w_i \\in \\{{\\phi_i\\}},\\quad {index_range(tests_used[0] + 1, tests_used[-1] + 1)}\\qquad \\text{{({len(tests_used)} denklem, {len(r["free"])} bilinmeyen)}}')
    for i in range(dofs):
        if i in tests_used:
            continue
        P = 'x_0' if i == 0 else 'x_1'
        if r['ivp']:
            D(f'w_{{{i + 1}}} = \\phi_{{{i + 1}}}\\ \\text{{çıkarıldı: }} {P}\\text{{ ucunda akı bilinmiyor, test fonksiyonu orada sıfır seçilir}}')
        else:
            D(f'w_{{{i + 1}}} = \\phi_{{{i + 1}}}\\ \\text{{çıkarıldı: }} u({P})\\ \\text{{verildi (temel koşul) }}\\Rightarrow a_{{{i + 1}}}\\ \\text{{bilinen; test uzayı: }} w({P}) = 0')

    def matrix_definitions():
        D(f'\\mathbf{{a}} = \\begin{{bmatrix}} a_1 & \\cdots & a_{{{dofs}}} \\end{{bmatrix}}^{{\\mathsf{{T}}}} \\in \\mathbb{{R}}^{{{dofs}}},\\quad '
          f'\\mathbf{{K}} \\in \\mathbb{{R}}^{{{dofs}\\times{dofs}}},\\quad \\mathbf{{F}} \\in \\mathbb{{R}}^{{{dofs}}}\\qquad '
          '\\text{(satır } i \\leftrightarrow w_i = \\phi_i\\text{, sütun } j \\leftrightarrow a_j\\text{; } K_{ij}\\text{: } a_j\\text{\'nin katsayısı)}')
        D('\\sum_{j} K_{ij}\\,a_j = F_i \\;\\Longleftrightarrow\\; \\mathbf{K}\\,\\mathbf{a} = \\mathbf{F}')

    if order == 1:
        D('\\text{Birinci mertebe: } R \\text{ yalnız } u_h\' \\text{ içerir, } C^0 \\text{ yaklaşımda integrallenebilir} \\Rightarrow \\text{kısmi integrasyon gerekmez}')
        for row in form['galerkin_rows']:
            parts = ' + '.join(f'\\int_{{{sp.latex(p_["xe"])}}}^{{{sp.latex(p_["xe"] + h)}}}\\left({sp.latex(p_["weight"])}\\right)\\left[{sp.latex(p_["residual"])}\\right]dx'
                               for p_ in row['parts'])
            D(f'w_{{{row["i"] + 1}}}:\\quad {parts} = {sp.latex(row["total"])} = 0')
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
                why = f"{name} = \\alpha({P})\\,u'({P}) = {tex_value(r['Q'][side])}\\ \\text{{(türev koşulu akıya çevrildi: α ile çarpıldı; yük vektörüne girer)}}"
            elif kind == 'flux':
                why = f"{name} = \\alpha({P})\\,u'({P}) = {tex_value(r['Q'][side])}\\ \\text{{(akı koşulu doğrudan; yük vektörüne girer)}}"
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
    D(f'K^{{(e)}}_{{ij}} = \\int_0^h\\left({integrand_K}\\right)ds,\\qquad F^{{(e)}}_{{l,i}} = {"-" if order == 2 else ""}\\int_0^h N_i\\,f(x_e + s)\\,ds,\\qquad x = x_e + s'
      + ('' if order == 2 else '\\qquad (\\mathbf{F} = \\mathbf{F}_l\\text{; sınır terimi yok})'))
    methods = {step['method'] for record in r['elements'] for step in record['K_steps'] + record['F_steps']}
    if not r['symbolic']:
        D('\\text{Sembolik süre sınırı aşıldı: tüm integraller sayısal (mpmath.quad, uyarlamalı tanh–sinh kuralı, 15 basamak); adımlar yazılmadı}')
    elif 'sayısal' in methods:
        D('\\approx:\\ \\text{kapalı biçim bulunamadı; sayısal integral (mpmath.quad, uyarlamalı tanh–sinh kuralı, 15 basamak)}')
    else:
        D('\\text{Gösterilen integraller SymPy ile tam (sembolik) alındı.}')
    if n > DETAIL_LIMIT and r['symbolic']:
        D(f'\\text{{gösterilmeyen elemanlar (e = 3, …, {n - 1}): sayısal integral (mpmath.quad, tanh–sinh); değerleri montajda}}')
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
            if r['symbolic']:
                for step in record['K_steps']:
                    D(format_integral(step))
            D(f'\\mathbf{{K}}^{{({e + 1})}} = {tex_matrix(record["K"].tolist())}')
        if has_F:
            if r['symbolic']:
                for step in record['F_steps']:
                    D(format_integral(step))
            D(f'\\mathbf{{F}}_l^{{({e + 1})}} = {tex_column(list(record["F"]))}')
    if order == 2 and not alpha.free_symbols and alpha < 0:
        first_K, first_F = r['elements'][0]['K'], r['elements'][0]['F']
        D(f'\\alpha = {sp.latex(alpha)} < 0:\\qquad \\tilde{{\\mathbf{{K}}}}^{{(1)}} = -\\mathbf{{K}}^{{(1)}} = {tex_matrix((-first_K).tolist())},'
          f'\\qquad \\tilde{{\\mathbf{{F}}}}_l^{{(1)}} = -\\mathbf{{F}}_l^{{(1)}} = {tex_column(list(-first_F))}')
        D('\\tilde{\\mathbf{K}}\\,\\mathbf{a} = \\tilde{\\mathbf{F}} \\iff (-\\mathbf{K})\\,\\mathbf{a} = -\\mathbf{F} \\iff '
          '\\mathbf{K}\\,\\mathbf{a} = \\mathbf{F}\\qquad \\text{(aynı }\\mathbf{a}\\text{; standart pozitif rijitlik işaretiyle)}')

    # --- montaj
    section('Montaj')
    D('K_{IJ} = \\sum_{e} K^{(e)}_{ij},\\qquad F_I = \\sum_e F^{(e)}_i,\\qquad I = (e-1)p + i,\\ J = (e-1)p + j\\qquad \\text{(yerel } i,j \\to \\text{ global } I,J\\text{)}')
    for c in form['contributions']:
        terms = ' + '.join(f'K^{{({e + 1})}}_{{{i + 1}{j + 1}}}' for e, i, j, _ in c['parts'])
        values = ' + '.join(f'\\left({tex(v)}\\right)' for *_, v in c['parts'])
        D(f'K_{{{c["I"] + 1}{c["J"] + 1}}} = {terms} = {values} = {tex(c["value"])}' + (f'\\qquad \\text{{(ortak düğüm {c["I"] + 1})}}' if c['I'] == c['J'] else ''))
    for c in form['load_contributions']:
        D(f'F_{{l,{c["I"] + 1}}} = ' + ' + '.join(f'F^{{({e + 1})}}_{{l,{i + 1}}}' for e, i, _ in c['parts']) + ' = '
          + ' + '.join(f'\\left({tex(v)}\\right)' for *_, v in c['parts']) + f' = {tex(c["value"])}')
    if form['contributions']:
        D('\\text{diğer girdiler tek bir elemandan gelir: } K_{IJ} = K^{(e)}_{ij},\\ \\text{boş girdiler: ortak elemanı olmayan düğümler}')
    D('\\mathbf{K} = ' + tex_matrix(r['K'].tolist()) + ',\\qquad \\mathbf{F}_l = ' + tex_column(list(r['F_l'])))
    if order == 2:
        D("\\mathbf{F}_b = \\begin{bmatrix} -Q_0 & 0 & \\cdots & 0 & Q_1 \\end{bmatrix}^{\\mathsf{T}}\\qquad \\text{(yalnız uç düğümler; iç düğümlerde } \\phi_i(x_0) = \\phi_i(x_1) = 0\\text{)}")
        for row in form['row_equations']:
            D(f'w_{{{row["i"] + 1}}}:\\quad {sp.latex(row["lhs"])} = {sp.latex(row["rhs"])}')

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
    if removed_rows:
        D('\\text{çıkarılan satırlar: } ' + ', '.join(str(i) for i in removed_rows) + '\\ (w_i \\text{ test uzayında değil})')
    if removed_cols:
        D('\\text{çıkarılan sütunlar: } ' + ', '.join(str(j) for j in removed_cols)
          + '\\ (a_j \\text{ bilinen} \\Rightarrow K_{ij}a_j \\text{ sağ tarafa: } F_{\\mathrm{red},i} = F_i - \\sum_{j\\ \\text{bilinen}} K_{ij}a_j)')
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
            flux = r['reactions'][side]
            slope = (f",\\qquad u_h'({sp.latex(point)}) = \\frac{{{name}}}{{\\alpha({sp.latex(point)})}} = {tex_value(r['reaction_slopes'][side])}"
                     if side in r['reaction_slopes'] else '')
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
    if form['ivp_slope']:
        sl = form['ivp_slope']
        D(f"u_h'({sp.latex(sl['point'])}) = {tex_value(sl['computed'])},\\qquad u'({sp.latex(sl['point'])}) = {sp.latex(sl['given'])}"
          f"\\qquad \\text{{(u' koşulu akı üzerinden, zayıf biçimde sağlanır)}}")

    # --- yeniden kurulum
    section('Yaklaşık çözümün yeniden kurulması')
    D('u_h^{(e)}(x) = \\sum_{k=1}^{p+1} \\hat{a}^{(e)}_k\\,N_k^{(e)}(x),\\qquad \\hat{a}^{(e)}_k = a_{(e-1)p+k}')
    for piece in form['reconstruction']:
        e, xe = piece['e'], piece['xe']
        if piece['N'] is not None:
            sub = ' + '.join(f'\\left({tex(sp.sympify(v))}\\right)\\left({sp.latex(Nk)}\\right)' for v, Nk in zip(piece['values'], piece['N']))
            D(f'u_h^{{({e + 1})}}(x) = {sub} = {tex(piece["piece"])},\\qquad {sp.latex(xe)} \\le x \\le {sp.latex(xe + h)}')
        else:
            D(f'u_h^{{({e + 1})}}(x) = {tex(piece["piece"])},\\qquad {sp.latex(xe)} \\le x \\le {sp.latex(xe + h)}')
    if n > DETAIL_LIMIT:
        D('\\text{diğer elemanlarda aynı biçimde}')

    # --- hata ölçüleri
    section('Hata ölçüleri')
    em = V['errors']
    if em is None:
        D('\\text{hesaplanmadı: doğrulanmış analitik çözüm yok (' + ('aday reddedildi' if exact['status'] == 'rejected' else exact['reason']) + ')}')
    else:
        D(f'u(x) = {sp.latex(solution)}\\qquad \\text{{({level[exact["status"]]})}}')
        D(',\\quad '.join(f'{name[1]}\\ \\checkmark\\ \\text{{({how})}}' for name, ok, how in exact['checks']))
        exact_or_dec = lambda v: sp.latex(v) if sp.sympify(v).is_Rational else dec(sp.sympify(v).evalf())
        table = ' \\\\ '.join(f'{j + 1} & {tex(xj)} & {exact_or_dec(aj)} & {exact_or_dec(uj)} & {exact_or_dec(ej)}'
                              for j, xj, aj, uj, ej in em['nodal_table'])
        D(f'\\begin{{array}}{{c|c|c|c|c}} j & x_j & a_j = u_h(x_j) & u(x_j) & |u - u_h| \\\\ \\hline {table} \\end{{array}}\\qquad \\text{{(düğüm hatası)}}')
        mids = ' \\\\ '.join(f'{tex(xm)} & {dec(um)} & {dec(uhm)} & {dec(um - uhm)}' for xm, um, uhm in em['midpoints'])
        D(f'\\begin{{array}}{{c|c|c|c}} x & u(x) & u_h(x) & u - u_h \\\\ \\hline {mids} \\end{{array}}\\qquad \\text{{(noktasal hata, her elemanda ilk iki düğümün ortası }} x_e + h/2p\\text{{)}}')
        rel = lambda v: dec(v) if v is not None else '\\text{tanımsız}'
        D(f'\\max_j |u(x_j) - a_j| = {dec(em["nodal"])}\\ (\\text{{göreli }} {rel(em["rel_nodal"])}),\\qquad '
          f'\\max_x |u - u_h| \\approx {dec(em["pointwise"])}\\ (\\text{{göreli }} {rel(em["rel_pointwise"])};\\ \\text{{eleman başına 200 örnek, alt tahmin}})')
        D(f'\\|u - u_h\\|_{{L^2}} = \\left(\\int_{{{X0}}}^{{{X1}}}(u - u_h)^2\\,dx\\right)^{{1/2}} \\approx {dec(em["l2"])}\\ (\\text{{göreli }} {rel(em["rel_l2"])}),\\qquad '
          f'|u - u_h|_{{H^1}} = \\left(\\int_{{{X0}}}^{{{X1}}}(u\' - u_h\')^2\\,dx\\right)^{{1/2}} \\approx {dec(em["h1"])}\\ (\\text{{göreli }} {rel(em["rel_h1"])})')
        if em['nodal'] < 1e-12 * max(1.0, abs(float(max(em['l2'], em['pointwise'])))) and em['pointwise'] > 1e-12:
            D('\\text{düğümlerde hata yok ama eleman içinde var: düğüm değerlerinin doğruluğu fonksiyonun doğruluğu değildir}')
    conv = V['convergence']
    if conv.get('rows'):
        table = ' \\\\ '.join(f"{row['n']} & {dec(row['h'])} & {dec(row['l2'])} & {dec(row['rate_l2']) if row.get('rate_l2') is not None else '-'} & "
                              f"{dec(row['h1'])} & {dec(row['rate_h1']) if row.get('rate_h1') is not None else '-'}" for row in conv['rows'])
        D(f'\\begin{{array}}{{c|c|c|c|c|c}} n & h & \\|e\\|_{{L^2}} & r_{{L^2}} & |e|_{{H^1}} & r_{{H^1}} \\\\ \\hline {table} \\end{{array}}'
          '\\qquad r = \\frac{\\log(E_h/E_{h/2})}{\\log 2}\\ \\text{(bağımsız çözücü ifem)}')
        expected = conv['expected']
        D('\\text{beklenen: }' + (f"L^2\\ {expected['l2']},\\ H^1\\ {expected['h1']}\\ \\text{{({expected['basis']})}}" if expected
                                 else '\\text{bu sınıf için kuramsal hız iddia edilmiyor}'))
        D('\\text{' + conv['detail'] + '}')
    else:
        D(f'\\text{{yakınsama: {conv["detail"]}}}')

    # --- doğrulama ve tanılar
    section('Doğrulama ve sayısal tanılar')
    for it in V['items']:
        D(f'{STATUS_MARK[it["status"]]}\\quad \\text{{{it["name"]} ({it["status"]}): {it["detail"]}}}')
    dist = V['distributional']
    if dist and dist.get('jumps') is not None and dist['jumps']:
        shown_j = dist['jumps'][:6]
        D(',\\ '.join(f'J_{{{(k + 1) * degree + 1}}} = {dec(v)}' for k, v in enumerate(shown_j)) + ('\\ \\dots' if len(dist['jumps']) > 6 else '')
          + '\\qquad \\text{(iç düğümlerde } \\alpha[u_h\'] \\text{; global numara)}')
        for side, (computed, given) in dist['ends'].items():
            P = X0 if side == 'left' else X1
            D(f"\\alpha u_h'({P}) = {dec(computed)},\\quad Q = {dec(given)}\\qquad \\text{{(doğal koşul ayrık çözümde yalnız zayıf biçimde sağlanır)}}")

    # --- uyarılar
    section('Uyarılar ve sınırlar')
    warnings = []
    g = r['degeneracy']
    if g['kind'] != 'uniform':
        warnings.append('önde gelen katsayı aralıkta sıfırlanıyor' + (' ve işaret değiştiriyor' if g['kind'] == 'sign_changing' else '')
                        + ': düzgün eliptiklik/koersivite yok; yakınsama kuramı ve kararlılık garanti değil')
    for it in V['items']:
        if it['status'] == 'uyarı':
            warnings.append(f"{it['name']}: {it['detail']}")
    if order == 1:
        warnings.append('birinci mertebe standart Galerkin kararlılık garantisi vermez; sert problemde salınım görülebilir')
    if r['ivp']:
        warnings.append("başlangıç değer probleminde u' koşulu zayıf biçimde sağlanır; kaba ağda salınımlı çözümde faz hatası birikir")
    if not r['symbolic'] or 'sayısal' in methods:
        warnings.append('bazı integraller sayısal alındı (≈); o değerler 15 basamaklıdır')
    unavailable = [it['name'] for it in V['items'] if it['status'] == 'yapılmadı']
    if unavailable:
        warnings.append('yapılmayan doğrulamalar: ' + ', '.join(unavailable))
    if V['errors'] is not None:
        warnings.append('noktasal en büyük hata örneklemeye dayanır (alt tahmin)')
    for w in warnings or ['özel uyarı yok']:
        D('\\text{' + w + '}')
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
        '\\begin{document}', '{\\LARGE\\textbf{Sonlu Eleman Çözüm Kâğıdı}}\\par\\smallskip\\hrule',
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
    return {'text': text, 'tex': build_tex(text), 'nodal': [float(sp.sympify(v).evalf()) for v in r['a']],
            'exact': str(solution) if solution is not None else None,
            'exact_js': exact_for_browser(solution, r['x0'], r['x1']), 'K': K, **verdict_payload(V['exact']),
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
        if error.level == 'sürekli problem':
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
        output, code = {'error': str(error)}, 1
    sys.stdout.write(json.dumps(output, ensure_ascii=False))
    sys.stdout.flush()
    os._exit(code)      # süre sınırını aşmış bir dsolve iş parçacığı varsa beklemeden çık


if __name__ == '__main__':
    import fem_sheet    # tek modül kopyası: fem_verify 'fem_sheet'i içe aktarır; __main__ ile ayrı kopya olmasın (Diagnosis aynı sınıf)
    fem_sheet.main()
