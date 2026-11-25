"""continuous_problem.py — SÜREKLİ problemin varlık / teklik / uyumluluk doğrulaması (dar kapsamlı, kanıt odaklı).

Neden ayrı: FEM Lab'ın ve scikit-fem'in sonuçları AYRIK sistem hakkındadır; tekil olmayan bir FEM matrisi sürekli problemin
çözümü olduğunu, tekil bir matris çözüm olmadığını kanıtlamaz. Bu modül uygulamanın kendi sınıflandırıcısını (fem_verify) ve
FEM motorunu ÇAĞIRMAZ; girdi metnini verification/ref_parser.expr ile okur (skfem_reference ile ortak) (üçüncü bir ayrıştırıcı yazılmaz).

Problem: (α u')' + β u' + γ u = f,  x ∈ [a, b];  2. mertebe ⇔ α ≢ 0.
  (R) düzenlilik: α, β, γ, f [a, b]'de sürekli (SymPy continuous_domain) ve önde gelen katsayının [a, b]'de sıfırı yok (solveset).
  Yöntem 1  başlangıç değer problemi (2. mertebe 'u ve u′'; 1. mertebe tek 'u ='): (R) altında doğrusal ODE varlık-teklik teoremi.
  Yöntem 2  düzenli iki noktalı SDP: doğrulanmış temel çözüm sistemi (L[y_i] = 0, L[u_p] = f sembolik; Wronskian ≠ 0) ve sınır
            matrisi B c = r. det B ≠ 0 → var ve tek; det B = 0: rank[B|r] > rank B → yok; aksi halde var ama tek değil.
  Yöntem 3  β ≡ 0 ve ayrık (separated) koşullar: Fredholm/Green uyumluluğu ∫ f v = [α v u' − α u v']ₐᵇ (çekirdek v için).
  Yöntem 4  doğrudan analitik kanıt: aday çözüm denkleme ve her koşula sembolik olarak konur (varlık; teklik değil).
Hiçbir zaman: çözücü başarısızlığından yokluk, sayısal örneklemeden varlık/teklik, FEM matrisinden sürekli hüküm.
"""
import threading

import mpmath
import sympy as sp

from ref_parser import expr, X

STATUSES = ('EXISTENCE_CONFIRMED', 'NONEXISTENCE_CONFIRMED', 'UNIQUENESS_CONFIRMED', 'NONUNIQUENESS_CONFIRMED',
            'COMPATIBILITY_SATISFIED', 'COMPATIBILITY_VIOLATED', 'INCONCLUSIVE', 'NOT_APPLICABLE')
ZERO_REL = 1e-30          # 50 basamakta |değer| ≤ 1e-30·ölçek → sıfır sayılmaz/sayılır kararı için eşik (evalf hata denetimli)
TIMEOUT = 25


def conclusion(status, method='', assumptions='', evidence=''):
    return {'status': status, 'method': method, 'assumptions': assumptions, 'evidence': evidence}


def is_zero(value, scale=1):
    """Sembolik sıfır kanıtı; olmazsa 50 basamaklı değerlendirme. Dönüş: (True/False/None, kanıt türü)."""
    try:
        if sp.simplify(value) == 0:
            return True, 'sembolik'
    except Exception:
        pass
    try:
        v = abs(complex(sp.N(value, 50)))
    except Exception:
        return None, 'değerlendirilemedi'
    s = max(1.0, abs(complex(sp.N(scale, 30))))
    if v > ZERO_REL * s:
        return False, f'50 basamak: |değer| = {v:.3e}'
    return None, f'50 basamak: |değer| = {v:.1e} (sıfır olduğu sembolik olarak gösterilemedi)'


def regularity(coeffs, lead, a, b):
    """(R)'yi kurar. Dönüş: (True/False/None, kanıt)."""
    interval = sp.Interval(a, b)
    for name, e in coeffs.items():
        try:
            dom = sp.calculus.util.continuous_domain(e, X, sp.S.Reals)
        except Exception:
            return None, f'{name}: süreklilik bölgesi belirlenemedi'
        if not interval.is_subset(dom):
            return False, f'{name} = {e}: [{a}, {b}] üzerinde sürekli değil'
    if not lead.free_symbols:
        return (lead != 0), ('önde gelen katsayı sabit ve sıfırdan farklı' if lead != 0 else 'önde gelen katsayı sıfır')
    try:
        zeros = sp.solveset(lead, X, interval)
    except Exception:
        return None, 'önde gelen katsayının sıfırları belirlenemedi'
    if zeros == sp.S.EmptySet:
        return True, f'önde gelen katsayı {lead} [{a}, {b}] üzerinde sıfırlanmıyor (solveset)'
    if isinstance(zeros, sp.FiniteSet):
        return False, f'önde gelen katsayı {lead} sıfırlanıyor: x ∈ {sorted(zeros, key=float)}'
    return None, f'önde gelen katsayının sıfırları kapalı biçimde bulunamadı ({type(zeros).__name__})'


def endpoint(e, xb, side):
    """Uç değer: tek taraflı limit (dejenere uçta x·log x gibi ifadeler için); sonlu değilse None."""
    try:
        v = sp.limit(e, X, xb, '+' if side == 'left' else '-')
    except Exception:
        return None
    return v if v.is_finite else None


def boundary_rows(cond, alpha, sides):
    """Her koşulu doğrusal fonksiyonel olarak verir: (ad, fonksiyonel(y) → değer, veri değeri)."""
    rows = []
    for side, xb in sides:
        kind, value, slope = cond[side]
        if kind in ('dirichlet', 'both'):
            rows.append((f'u({xb})', lambda y, xb=xb, side=side: endpoint(y, xb, side), expr(value)))
        if kind == 'neumann':
            rows.append((f"u'({xb})", lambda y, xb=xb, side=side: endpoint(sp.diff(y, X), xb, side), expr(value)))
        if kind == 'both':
            rows.append((f"u'({xb})", lambda y, xb=xb, side=side: endpoint(sp.diff(y, X), xb, side), expr(slope)))
        if kind == 'flux':
            rows.append((f"αu'({xb})", lambda y, xb=xb, side=side: endpoint(alpha * sp.diff(y, X), xb, side), expr(value)))
        if kind == 'none':
            rows.append((f"αu'({xb}) (doğal)", lambda y, xb=xb, side=side: endpoint(alpha * sp.diff(y, X), xb, side), sp.Integer(0)))
    return rows


def fundamental_system(alpha, beta, gamma, f, a, b):
    """dsolve → (y1, y2, u_p) ve bağımsız doğrulama. Başarısızlık INCONCLUSIVE'dir (yokluk değil)."""
    u = sp.Function('u')
    L = lambda y: sp.diff(alpha * sp.diff(y, X), X) + beta * sp.diff(y, X) + gamma * y
    try:
        general = sp.dsolve(sp.Eq(L(u(X)), f), u(X)).rhs
    except NotImplementedError:
        return None, 'SymPy dsolve kapalı biçimde çözüm bulamadı (bu, çözümün olmadığı anlamına gelmez)'
    constants = sorted(general.free_symbols - {X}, key=str)
    if general.has(sp.Integral) or len(constants) != 2:
        return None, 'kapalı biçimde iki sabitli genel çözüm bulunamadı'
    ys = [sp.diff(general, c) for c in constants]
    up = general.subs({c: 0 for c in constants})
    if any(y.free_symbols - {X} for y in ys) or sp.simplify(general - up - sum(c * y for c, y in zip(constants, ys))) != 0:
        return None, 'genel çözüm sabitlerde doğrusal değil'
    for name, residual in (('L[y1]', L(ys[0])), ('L[y2]', L(ys[1])), ('L[u_p] − f', L(up) - f)):
        if sp.simplify(residual) != 0:
            return None, f'{name} = 0 sembolik olarak doğrulanamadı'
    # Doğrusal bağımsızlık: W = y1 y2′ − y2 y1′ bir x* noktasında ≠ 0 (Abel: α W = sabit, α ≠ 0 olan aralıkta hiç sıfırlanmaz).
    # Bağımlı bir çift det B ≡ 0 verir ve sahte "tek değil / çözüm yok" hükmüne yol açardı.
    W = ys[0] * sp.diff(ys[1], X) - ys[1] * sp.diff(ys[0], X)
    for t in (sp.Rational(1, 2), sp.Rational(1, 3), sp.Rational(2, 3)):
        xs = a + (b - a) * t
        try:
            w = W.subs(X, xs)
        except Exception:
            continue
        nonzero, how = is_zero(w, 1)
        if nonzero is False:
            return (ys, up, L), (f'L[y1] = L[y2] = 0 ve L[u_p] = f sembolik olarak doğrulandı; W({sp.sstr(xs)}) ≠ 0 ({how})')
    return None, 'temel çözümlerin bağımsızlığı (Wronskian ≠ 0) gösterilemedi'


# x gerçekken arg(g + i·c) = atan2(c, g) ve sabit c ≠ 0 için atan2(c, g) = sign(c)·π/2 − atan(g/c) (her gerçek g için tam).
def real_rewrite(e):
    rule = lambda y, g: sp.sign(y) * sp.pi / 2 - sp.atan(g / y) if (y.is_number and y != 0) else sp.atan2(y, g)
    return e.replace(sp.arg, lambda z: rule(sp.im(z), sp.re(z))).replace(sp.atan2, rule)


def numeric_rank(M):
    if M.rows == 0:
        return 0
    with mpmath.workdps(60):
        A = mpmath.matrix(M.evalf(60).tolist())
        sv = mpmath.svd_c(A, compute_uv=False)
        top = max([abs(v) for v in sv] + [mpmath.mpf(1)])
        return sum(1 for v in sv if abs(v) > mpmath.mpf('1e-40') * top)


def analyze(inp):
    box = {}

    def work():
        try:
            box['res'] = _analyze(inp)
        except Exception as error:                     # beklenmeyen hata: kanıt yok → INCONCLUSIVE, nedeni yazılı
            box['res'] = {'class': '?', 'existence': conclusion('INCONCLUSIVE', evidence=f'{type(error).__name__}: {error}'),
                          'uniqueness': conclusion('INCONCLUSIVE'), 'compatibility': conclusion('INCONCLUSIVE')}
    thread = threading.Thread(target=work, daemon=True)
    thread.start()
    thread.join(TIMEOUT)
    return box.get('res', {'class': '?', 'existence': conclusion('INCONCLUSIVE', evidence=f'zaman aşımı ({TIMEOUT} s)'),
                           'uniqueness': conclusion('INCONCLUSIVE', evidence='zaman aşımı'),
                           'compatibility': conclusion('INCONCLUSIVE', evidence='zaman aşımı')})


def _analyze(inp):
    alpha, beta, gamma, f = (expr(inp.get(k)) for k in ('alpha', 'beta', 'gamma', 'f'))
    a, b = expr(inp['x0']), expr(inp['x1'])
    order = 2 if alpha != 0 else 1
    cond = {s: (inp[s].get('type', 'none'), inp[s].get('value', ''), inp[s].get('slope', '')) for s in ('left', 'right')}
    ivp = order == 1 or any(cond[s][0] == 'both' for s in cond)
    lead = alpha if order == 2 else beta
    reg, reg_evidence = regularity({'α': alpha, 'β': beta, 'γ': gamma, 'f': f}, lead, a, b)
    out = {'class': ('IVP' if ivp else 'BVP') + f', mertebe {order}', 'regularity': {'established': reg, 'evidence': reg_evidence}}
    NA = conclusion('NOT_APPLICABLE')

    if ivp:                                              # Yöntem 1
        # Yöntem 1 yalnız TAM bir başlangıç koşulu takımı için geçerlidir: 1. mertebe → bir uçta 'u =', öbür uç verilmemiş;
        # 2. mertebe → bir uçta 'u ve u′', öbür uç verilmemiş. Fazlası aşırı belirli (genelde çözümsüz), eksiği belirsiz.
        kinds = sorted(cond[s][0] for s in cond)
        scope = kinds == (['dirichlet', 'none'] if order == 1 else ['both', 'none'])
        if not scope:
            why = f'koşul takımı {kinds} Yöntem 1 kapsamında değil (tam bir başlangıç koşulu takımı gerekir)'
            out.update(existence=conclusion('INCONCLUSIVE', evidence=why), uniqueness=conclusion('INCONCLUSIVE', evidence=why),
                       compatibility=conclusion('NOT_APPLICABLE'))
            return out
        if reg:
            thm = 'doğrusal ODE için varlık-teklik teoremi (Picard–Lindelöf, global; doğrusal sistem)'
            out['existence'] = conclusion('EXISTENCE_CONFIRMED', 'Yöntem 1: ' + thm, '(R): ' + reg_evidence, 'tek koşul noktasında başlangıç değeri verilmiş')
            out['uniqueness'] = conclusion('UNIQUENESS_CONFIRMED', 'Yöntem 1: ' + thm, '(R): ' + reg_evidence, '')
        else:
            out['existence'] = conclusion('INCONCLUSIVE', evidence='(R) kurulamadı: ' + reg_evidence)
            out['uniqueness'] = conclusion('INCONCLUSIVE', evidence='(R) kurulamadı: ' + reg_evidence)
        out['compatibility'] = NA
        return out

    sides = (('left', a), ('right', b))
    fs, why = fundamental_system(alpha, beta, gamma, f, a, b)
    if fs is None:
        reason = f'Yöntem 2 uygulanamadı: {why}'
        out.update(existence=conclusion('INCONCLUSIVE', evidence=reason), uniqueness=conclusion('INCONCLUSIVE', evidence=reason),
                   compatibility=conclusion('INCONCLUSIVE', evidence=reason))
        return out
    ys, up, L = fs
    rows = boundary_rows(cond, alpha, sides)
    B, r = [], []
    for name, functional, data in rows:
        entries = [functional(y) for y in ys]
        u0 = functional(up)
        if any(v is None for v in entries + [u0]):
            reason = f'{name}: temel çözümlerden biri bu uçta sonlu değil (tekil uç); sınır matrisi kurulamadı'
            out.update(existence=conclusion('INCONCLUSIVE', evidence=reason), uniqueness=conclusion('INCONCLUSIVE', evidence=reason),
                       compatibility=conclusion('INCONCLUSIVE', evidence=reason))
            return out
        B.append(entries)
        r.append(data - u0)
    B, r = sp.Matrix(B), sp.Matrix(r)
    det = sp.simplify(B.det())
    scale = max([1] + [abs(complex(sp.N(v, 20))) for v in B]) ** 2
    zero, how = is_zero(det, scale)
    evidence_B = f'{why}; B = {sp.sstr(B.tolist())}, r = {sp.sstr(list(r))}, det B = {sp.sstr(det)} ({how})'
    method2 = 'Yöntem 2: doğrulanmış temel çözüm sistemi + sınır matrisi ' + ('(düzenli SDP)' if reg else '(R kurulamadı: yalnız doğrudan kanıt)')

    # Yöntem 4 (doğrudan kanıt): det B ≠ 0 ise tek aday; sembolik olarak denkleme ve koşullara konur
    def direct_candidate(c):
        cand = sp.simplify(up + sum(ci * y for ci, y in zip(c, ys)))
        # Klasik gerçek çözüm: [a, b]'nin tamamında sonlu ve sürekli, gerçek değerli (Im u ≡ 0); aksi halde kanıt sayılmaz.
        if cand.has(sp.I):                         # karmaşık biçim: önce Im u ≡ 0 (x gerçek) gösterilir, sonra gerçek kısım
            im_part = sp.simplify(real_rewrite(sp.im(sp.expand_complex(cand))))
            if im_part != 0:
                return None, f'aday u = {sp.sstr(cand)} gerçek değerli gösterilemedi (Im u = {sp.sstr(im_part)})'
            cand = sp.simplify(real_rewrite(sp.re(sp.expand_complex(cand))))
        try:
            dom = sp.calculus.util.continuous_domain(cand, X, sp.S.Reals)
        except Exception:
            return None, 'adayın süreklilik bölgesi belirlenemedi'
        if not sp.Interval(a, b).is_subset(dom):
            return None, f'aday u = {sp.sstr(cand)} [{a}, {b}] üzerinde sürekli/sonlu değil'
        if sp.simplify(L(cand) - f) != 0:
            return None, 'aday denklemi sembolik olarak sağlamadı'
        for (name, functional, data) in rows:
            v = functional(cand)
            if v is None or is_zero(v - data, 1)[0] is not True:
                return None, f'aday {name} koşulunu sağlamadı'
        return cand, 'aday denkleme ve her koşula sembolik olarak kondu'

    if zero is False:
        c = list(B.LUsolve(r))
        cand, cand_ev = direct_candidate(c)
        if reg:
            out['existence'] = conclusion('EXISTENCE_CONFIRMED', method2, '(R): ' + reg_evidence, evidence_B + '; ' + cand_ev)
            out['uniqueness'] = conclusion('UNIQUENESS_CONFIRMED', method2, '(R): ' + reg_evidence,
                                           'det B ≠ 0: homojen problemin yalnız sıfır çözümü var')
        else:
            out['existence'] = (conclusion('EXISTENCE_CONFIRMED', 'Yöntem 4: doğrudan analitik kanıt', 'düzenlilik gerekmez (çözüm açıkça sınandı)',
                                           f'u = {sp.sstr(cand)}; ' + cand_ev) if cand is not None
                                else conclusion('INCONCLUSIVE', evidence='(R) yok ve ' + cand_ev))
            out['uniqueness'] = conclusion('INCONCLUSIVE', evidence='(R) kurulamadı (' + reg_evidence + '): tekil problemde teklik kuramı uygulanmaz')
        out['compatibility'] = (conclusion('NOT_APPLICABLE', evidence='det B ≠ 0: çekirdek yalnız sıfır, uyumluluk koşulu yok') if reg
                                else conclusion('INCONCLUSIVE', evidence='(R) yok: temel sistem aralığın tamamında geçerli olmayabilir; '
                                                'det B ≠ 0 çekirdeğin sıfır olduğunu kanıtlamaz'))
        out['solution'] = sp.sstr(cand) if cand is not None else None
        return out
    if zero is None or not reg:
        reason = ('det B = 0 kanıtlanamadı ve ≠ 0 da gösterilemedi: ' + how) if zero is None else '(R) kurulamadı: ' + reg_evidence
        out.update(existence=conclusion('INCONCLUSIVE', evidence=reason), uniqueness=conclusion('INCONCLUSIVE', evidence=reason),
                   compatibility=conclusion('INCONCLUSIVE', evidence=reason))
        return out

    # det B = 0 (sembolik) ve (R): çekirdek sıfırdan farklı
    rank_B, rank_aug = B.rank(simplify=True), B.row_join(r).rank(simplify=True)
    num = (numeric_rank(B), numeric_rank(B.row_join(r)))
    rank_ev = f'rank B = {rank_B}, rank [B|r] = {rank_aug} (sembolik; 60 basamaklı karmaşık SVD: {num[0]}, {num[1]})'
    if (rank_B, rank_aug) != num:
        out.update(existence=conclusion('INCONCLUSIVE', evidence='sembolik ve sayısal rank uyuşmuyor: ' + rank_ev),
                   uniqueness=conclusion('INCONCLUSIVE', evidence=rank_ev), compatibility=conclusion('INCONCLUSIVE', evidence=rank_ev))
        return out
    kernel = [sp.simplify(sum(vi * y for vi, y in zip(v, ys))) for v in B.nullspace()]
    if rank_aug > rank_B:
        out['existence'] = conclusion('NONEXISTENCE_CONFIRMED', method2, '(R): ' + reg_evidence,
                                      evidence_B + '; ' + rank_ev + ': r, B\'nin değer uzayında değil → hiçbir u_p + c·y koşulları sağlamaz')
        out['uniqueness'] = conclusion('NOT_APPLICABLE', evidence='çözüm yok')
    else:
        out['existence'] = conclusion('EXISTENCE_CONFIRMED', method2, '(R): ' + reg_evidence, evidence_B + '; ' + rank_ev)
        out['uniqueness'] = conclusion('NONUNIQUENESS_CONFIRMED', method2, '(R): ' + reg_evidence,
                                       f'homojen problemin sıfırdan farklı çözümü: v = {", ".join(sp.sstr(v) for v in kernel)}')

    # Yöntem 3: β ≡ 0 ve ayrık koşullar → bağımsız Fredholm/Green uyumluluğu (Yöntem 2 ile çapraz kontrol)
    if beta != 0:
        out['compatibility'] = conclusion('NOT_APPLICABLE', evidence='β ≢ 0: operatör biçimsel olarak öz-eşlenik değil; Yöntem 3 kullanılmaz')
        return out
    results = []
    for v in kernel:
        terms = []
        for side, xb in sides:
            kind, value, _ = cond[side]
            if kind == 'dirichlet':              # v(x_b) = 0, bilinmeyen u' düşer
                term = -endpoint(alpha * sp.diff(v, X), xb, side) * expr(value)
            elif kind == 'neumann':              # α v' = 0, bilinmeyen u düşer
                term = endpoint(alpha * v, xb, side) * expr(value)
            elif kind == 'flux':
                term = endpoint(v, xb, side) * expr(value)
            else:                                # doğal: α u' = 0 ve α v' = 0
                term = sp.Integer(0)
            terms.append(term)
        bracket = terms[1] - terms[0]
        integral = sp.integrate(sp.expand(f * v), (X, a, b))
        if integral.has(sp.Integral):
            with mpmath.workdps(50):
                integral = sp.Float(mpmath.quad(sp.lambdify(X, f * v, 'mpmath'), [float(a), float(b)]), 40)
            how_int = 'sayısal kuadratür (50 basamak)'
        else:
            how_int = 'sembolik integral'
        C = sp.simplify(integral - bracket)
        z, how_c = is_zero(C, 1)
        results.append((v, integral, bracket, C, z, how_int, how_c))
    if any(z is None for *_, z, _, _ in results):
        out['compatibility'] = conclusion('INCONCLUSIVE', 'Yöntem 3: Fredholm/Green', 'β ≡ 0, ayrık koşullar',
                                          '; '.join(f'C = {sp.sstr(C)} ({hc})' for _, _, _, C, _, _, hc in results))
        return out
    violated = any(z is False for *_, z, _, _ in results)
    detail = '; '.join(f'v = {sp.sstr(v)}: ∫ f v dx = {sp.sstr(I)} ({hi}), [α v u′ − α u v′]ₐᵇ = {sp.sstr(Bk)}, fark C = {sp.sstr(C)} ({hc})'
                       for v, I, Bk, C, z, hi, hc in results)
    out['compatibility'] = conclusion('COMPATIBILITY_VIOLATED' if violated else 'COMPATIBILITY_SATISFIED',
                                      'Yöntem 3: Fredholm alternatifi, Green özdeşliği ∫ v L u − ∫ u L v = [α(v u′ − u v′)]ₐᵇ',
                                      'β ≡ 0 (L biçimsel olarak öz-eşlenik), ayrık sınır koşulları, (R)', detail)
    agree = violated == (out['existence']['status'] == 'NONEXISTENCE_CONFIRMED')
    out['compatibility']['cross_check'] = 'Yöntem 2 ile tutarlı' if agree else 'Yöntem 2 ile ÇELİŞKİLİ'
    return out
