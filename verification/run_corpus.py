#!/usr/bin/env python3
"""run_corpus.py — 80 vakalık YALNIZ GİRDİ korpusunu gerçek uygulama hattından geçirir ve sonuçları bağımsız olarak değerlendirir.

    python3 verification/run_corpus.py                    (tümü)
    python3 verification/run_corpus.py C12 C38            (seçili vakalar)

Korpus (verification/corpus/cases.json) beklenen sonuç İÇERMEZ; her vaka {id, meta, input}. Motora yalnız 'input' gider.
Sonuçlar verification/corpus_results/ altına yazılır (korpustan ayrı): results.json, skfem_results.json, report.md.

Her vaka için kullanılan yollar (hiçbiri atlanmaz, doğrulama bypass edilmez):
  1. Sunucu hattı: server.compute_sheet(input) — tarayıcının POST /api/sheet ile çağırdığı yolun aynısı (ayrı süreç, süre
     sınırı, sembolik → sayısal yedek). 200 = kâğıt, 400 = temiz ret, 500 = çöküş, 504 = zaman aşımı.
  2. Tarayıcı çekirdeği: node tests.js --solve (index.html CORE: buildProblem + METHODS.galerkin.solve).
  3. scikit-fem bağımsız referansı (verification/skfem_reference.py) — zayıf formu kendisi montajlar.
  4. Analitik çözüm, kâğıt bulduysa, BURADA ayrıca (bağımsız ayrıştırıcıyla) denklem ve koşullara karşı sınanır; hata
     normları scikit-fem dışı, Gauss–Legendre ile yeniden hesaplanıp kâğıtta yazanla karşılaştırılır.

Durumlar (vaka): PASS · WARN (sayısal uyarı / kuadratürle sınırlı) · LIMITATION (belgelenmiş sınır) · REJECTED_OK (geçersiz
girdi temiz reddedildi) · DIAGNOSED (geçerli ama tekil/uyumsuz: teşhis bağımsız kanıtla tutarlı) · FAIL · UNAVAILABLE.
Toleranslar tek bir evrensel sayı değildir; her karşılaştırmanın yanında gerekçesiyle yazılır.
"""
import json
import os
import re
import subprocess
import sys
import time
from pathlib import Path

import numpy as np
import sympy as sp

HERE = Path(__file__).resolve().parent
APP = HERE.parent
sys.path.insert(0, str(APP))
sys.path.insert(0, str(HERE))
import server                    # noqa: E402  gerçek sunucu hattı (compute_sheet)
import skfem_reference as SK     # noqa: E402  bağımsız referans
import continuous_problem as CP  # noqa: E402  sürekli problem doğrulaması (uygulamadan bağımsız)

EPS = np.finfo(float).eps
CORPUS = HERE / 'corpus' / 'cases.json'
OUT = HERE / 'corpus_results'
JS_CRASH = ('TypeError', 'RangeError', 'ReferenceError', 'SyntaxError', 'InternalError')


def js_solve(inp):
    """Tarayıcı çekirdeği. Dönüş: ('ok', c) · ('rejected', mesaj) · ('crash', mesaj)."""
    run = subprocess.run(['node', str(APP / 'tests.js'), '--solve', json.dumps(inp)], capture_output=True, text=True, timeout=120)
    if run.returncode == 0:
        return 'ok', json.loads(run.stdout)
    # node yakalanmamış hatayı "<Tür>: mesaj" satırıyla yazar (Diagnosis için "Diagnosis [Error]: …"); kaynak satırı değil o alınır
    line = next((l for l in run.stderr.splitlines() if re.match(r'^\w*(Error|Diagnosis)( \[\w+\])?: ', l)), run.stderr.strip()[-200:])
    kind = line.split(':')[0].strip()
    return ('crash' if any(k in kind for k in JS_CRASH) else 'rejected'), line[:240]


def scale_tol(base, cond):
    """Göreli tolerans: temel değer ile koşul sayısının 64-bit yuvarlamaya etkisinin büyüğü (κ·ε, 100 kat payla)."""
    return max(base, 100.0 * (cond or 1.0) * EPS)


def rel_gap(a, b):
    a, b = np.asarray(a, dtype=float), np.asarray(b, dtype=float)
    return float(np.max(np.abs(a - b)) / max(1.0, float(np.max(np.abs(b)))))


# Kâğıdın bulduğu analitik çözümü bağımsız ayrıştırıcıyla okuyup denklem ve koşullara karşı sınar (uygulamanın verify_exact'i
# kullanılmaz). Dönüş: (geçerli mi, kanıt).
def check_exact(inp, text):
    try:
        u = SK.expr(text.replace('**', '^'))
        al, be, ga, f = (SK.expr(inp.get(k)) for k in ('alpha', 'beta', 'gamma', 'f'))
    except Exception as error:
        return None, f'okunamadı ({type(error).__name__})'
    x = SK.X
    R = sp.diff(al * sp.diff(u, x), x) + be * sp.diff(u, x) + ga * u - f
    x0, x1 = float(SK.expr(inp['x0'])), float(SK.expr(inp['x1']))
    pts = np.linspace(x0, x1, 23)[1:-1]
    Rf = sp.lambdify(x, R, 'numpy')
    scale = max(1.0, float(np.max(np.abs(sp.lambdify(x, u, 'numpy')(pts) * np.ones_like(pts)))))
    res = float(np.max(np.abs(Rf(pts) * np.ones_like(pts))))
    bad = []
    if not res <= 1e-8 * scale:
        bad.append(f'|R[u]| = {res:.1e}')
    du = sp.diff(u, x)
    for s, xb in (('left', x0), ('right', x1)):
        c = inp[s]
        kind = c.get('type')
        val = lambda t: float(SK.expr(t))
        ub, dub = float(u.subs(x, xb)), float(du.subs(x, xb))
        if kind in ('dirichlet', 'both') and abs(ub - val(c['value'])) > 1e-9 * scale:
            bad.append(f'u({xb}) = {ub} ≠ {c["value"]}')
        if kind == 'neumann' and abs(dub - val(c['value'])) > 1e-9 * scale:
            bad.append(f"u'({xb}) = {dub} ≠ {c['value']}")
        if kind == 'both' and abs(dub - val(c['slope'])) > 1e-9 * scale:
            bad.append(f"u'({xb}) = {dub} ≠ {c['slope']}")
        if kind == 'flux' and abs(float(al.subs(x, xb)) * dub - val(c['value'])) > 1e-9 * scale:
            bad.append(f"α u'({xb}) ≠ {c['value']}")
    return (not bad), ('denklem ve koşullar bağımsız olarak sağlandı' if not bad else '; '.join(bad))


# Hata normları, uygulamadan bağımsız: FEM Lab'ın düğüm katsayılarından Lagrange yeniden kurulumu + eleman başına 40 noktalı
# Gauss–Legendre (yüksek dereceli düzgün çözüm için yeterli; kuadratür hatası ≪ 1e-10 göreli).
def error_norms(inp, a, exact_text):
    u = SK.expr(exact_text.replace('**', '^'))
    uf, duf = sp.lambdify(SK.X, u, 'numpy'), sp.lambdify(SK.X, sp.diff(u, SK.X), 'numpy')
    x0, x1, n, p = float(SK.expr(inp['x0'])), float(SK.expr(inp['x1'])), int(inp['n']), int(inp['degree'])
    t, w = np.polynomial.legendre.leggauss(40)
    h = (x1 - x0) / n
    l2 = h1 = 0.0
    for e in range(n):
        xs = x0 + e * h + h * (t + 1) / 2
        uh, duh = SK.evaluate(a, x0, x1, n, p, np.clip(xs, x0, x1 - 1e-15 * max(1, abs(x1))))
        l2 += float(np.sum(w * h / 2 * (uf(xs) * np.ones_like(xs) - uh) ** 2))
        h1 += float(np.sum(w * h / 2 * (duf(xs) * np.ones_like(xs) - duh) ** 2))
    return l2 ** 0.5, h1 ** 0.5


def sheet_number(text, pattern):
    m = re.search(pattern, text)
    if not m:
        return None
    raw = m.group(1)
    m2 = re.match(r'(-?[\d.]+)\\times 10\^\{(-?\d+)\}', raw)
    return float(m2.group(1)) * 10 ** int(m2.group(2)) if m2 else float(raw)


MATRIX = ('tarayıcı ↔ sunucu', 'FEM Lab ↔ ifem', 'scikit-fem eleman', 'scikit-fem global K', 'scikit-fem düğüm', 'analitik çözüm',
          'sürekli varlık', 'sürekli teklik', 'uyumluluk', 'L²/H¹ yeniden hesap', 'sınır koşulu')
APP_ITEM = {'passed': 'PASS', 'warning': 'WARN', 'failed': 'FAIL', 'not_run': 'NOT_RUN'}   # uygulamanın makine durum değerleri


def app_verdict(code, sheet):
    """Uygulamanın sürekli düzey hükmü: 'solved' · 'K1' · 'K2' (sürekli) · 'discrete' (yalnız ayrık hüküm) · 'other'."""
    if code == 200:
        return 'solved'
    sheet = sheet if isinstance(sheet, dict) else {}
    # Uygulamanın makine alanları (fem_sheet.main: category, level). Mesaj metni ayrıştırılmaz: arayüz dili değişse de
    # hüküm aynı okunur (önceki sürüm "— ayrık sistem:" gibi Türkçe işaretlere bakıyordu).
    level, category = sheet.get('level'), sheet.get('category')
    if level == 'discrete':
        return 'discrete'
    if level == 'continuous':
        return category if category in ('K1', 'K2') else 'other'
    return 'other'


# Sürekli düzey sonucu ile uygulamanın davranışının tutarlılığı. Kesin (CONFIRMED) bir sürekli sonuç uygulamanın davranışıyla
# çelişirse FAIL; INCONCLUSIVE hiçbir zaman PASS sayılmaz.
def continuous_cells(cont, verdict):
    ex, un, co = cont['existence']['status'], cont['uniqueness']['status'], cont['compatibility']['status']
    def cell(status, ok_when, bad_when):
        if status in ('INCONCLUSIVE',):
            return 'INCONCLUSIVE'
        if status in ('NOT_APPLICABLE',):
            return 'NOT_RUN'
        if status in bad_when:
            return 'FAIL'
        return 'PASS' if status in ok_when else 'INCONCLUSIVE'
    if verdict == 'solved':
        e = cell(ex, ('EXISTENCE_CONFIRMED',), ('NONEXISTENCE_CONFIRMED',))
        u = cell(un, ('UNIQUENESS_CONFIRMED',), ('NONUNIQUENESS_CONFIRMED',))
        c = cell(co, ('COMPATIBILITY_SATISFIED',), ('COMPATIBILITY_VIOLATED',))
    elif verdict == 'K1':
        e = cell(ex, ('NONEXISTENCE_CONFIRMED',), ('EXISTENCE_CONFIRMED',))
        u = 'NOT_RUN' if un == 'NOT_APPLICABLE' else cell(un, (), ())
        c = cell(co, ('COMPATIBILITY_VIOLATED',), ('COMPATIBILITY_SATISFIED',))
    elif verdict == 'K2':
        e = cell(ex, ('EXISTENCE_CONFIRMED',), ('NONEXISTENCE_CONFIRMED',))
        u = cell(un, ('NONUNIQUENESS_CONFIRMED',), ('UNIQUENESS_CONFIRMED',))
        c = cell(co, ('COMPATIBILITY_SATISFIED',), ('COMPATIBILITY_VIOLATED',))
    else:                                              # ayrık hüküm: sürekli sonuç bilgi amaçlı, çelişki kavramı yok
        e = 'INCONCLUSIVE' if ex == 'INCONCLUSIVE' else 'PASS'
        u = 'INCONCLUSIVE' if un == 'INCONCLUSIVE' else 'NOT_RUN' if un == 'NOT_APPLICABLE' else 'PASS'
        c = 'INCONCLUSIVE' if co == 'INCONCLUSIVE' else 'NOT_RUN' if co == 'NOT_APPLICABLE' else 'PASS'
    return e, u, c


def element_comparison(sheet_K, sheet_F, ref, qtol):
    """Eleman eleman: FEM Lab'ın montajladığı K^(e), F^(e) ↔ scikit-fem tek eleman montajı (Φ ile düğümsel tabana çevrilmiş)."""
    worst = {'rel': 0.0, 'abs': 0.0, 'where': None}
    for kind, mine, theirs in (('K', sheet_K, ref['element_K']), ('F', sheet_F, ref['element_F'])):
        for e, (A, B) in enumerate(zip(mine, theirs)):
            A, B = np.asarray(A, dtype=float), np.asarray(B, dtype=float)
            diff = np.abs(A - B)
            scale = max(1.0, float(np.max(np.abs(A))))
            idx = np.unravel_index(int(np.argmax(diff)), diff.shape)
            rel = float(diff[idx] / scale)
            if rel > worst['rel']:
                worst = {'rel': rel, 'abs': float(diff[idx]), 'where': f'{kind}^({e + 1}) girdi {tuple(int(i) + 1 for i in idx)}'}
    ok = worst['rel'] <= qtol
    return ('PASS' if ok else 'FAIL'), (f"{len(mine)} eleman; en büyük göreli fark {worst['rel']:.1e} (mutlak {worst['abs']:.1e})"
                                        + (f" @ {worst['where']}" if worst['where'] else '') + f" (tolerans {qtol:.0e})"), worst


def run_case(case):
    inp, meta = case['input'], case['meta']
    rec = {'id': case['id'], 'name': meta['name'], 'category': meta['category'], 'validity': meta['validity'],
           'checks': [], 'matrix': {k: 'NOT_RUN' for k in MATRIX}}
    add = lambda name, status, evidence, layer=None: rec['checks'].append(
        {'check': name, 'status': status, 'evidence': evidence, **({'layer': layer} if layer else {})})
    M = rec['matrix']
    t0 = time.time()
    try:
        sheet, code = server.compute_sheet(inp)
    except Exception as error:
        sheet, code = {'error': f'{type(error).__name__}: {error}'}, 'exception'
    rec['server'] = {'code': code, 'seconds': round(time.time() - t0, 1), 'error': sheet.get('error')}
    try:
        js_kind, js_val = js_solve(inp)
    except subprocess.TimeoutExpired:
        js_kind, js_val = 'crash', 'zaman aşımı (120 s)'
    rec['js'] = {'kind': js_kind, 'detail': js_val if js_kind != 'ok' else None}

    if meta['validity'] == 'invalid':
        rec['matrix'] = {k: 'NOT_RUN' for k in MATRIX}
        ok_py = code == 400 and bool(sheet.get('error'))
        add('sunucu temiz ret', 'PASS' if ok_py else 'FAIL', f'HTTP {code}: {sheet.get("error") or "(kâğıt üretildi)"}',
            None if ok_py else 'input validation (sunucu)')
        if js_kind == 'rejected':
            add('tarayıcı çekirdeği temiz ret', 'PASS', js_val)
        else:
            add('tarayıcı çekirdeği temiz ret', 'FAIL', f'{js_kind}: {js_val if js_kind != "ok" else "girdi kabul edildi, c = " + json.dumps(js_val)[:60]}',
                'input validation (tarayıcı çekirdeği)')
        rec['status'] = 'REJECTED_OK' if all(c['status'] == 'PASS' for c in rec['checks']) else 'FAIL'
        return rec, None, None

    ref = SK.reference(inp)
    skrec = {'id': case['id'], 'status': ref['status'], 'skfem_version': ref['skfem_version'],
             **{k: ref.get(k) for k in ('element', 'degree', 'n', 'intorder', 'intorder_check', 'quadrature_estimate',
                                        'formulation', 'cond', 'solver_residual', 'reason')}}
    cont = CP.analyze(inp)
    contrec = {'id': case['id'], **cont}
    verdict = app_verdict(code, sheet)
    M['sürekli varlık'], M['sürekli teklik'], M['uyumluluk'] = continuous_cells(cont, verdict)
    for key, part in (('sürekli varlık', 'existence'), ('sürekli teklik', 'uniqueness'), ('uyumluluk', 'compatibility')):
        c = cont[part]
        add(f'sürekli: {part}', M[key], f"{c['status']} — {c.get('method', '')} {c.get('evidence', '')}"[:400],
            'continuous verification ↔ uygulama hükmü' if M[key] == 'FAIL' else None)

    if code in (500, 504, 'exception'):
        add('sunucu hattı', 'FAIL', f'HTTP {code}: {sheet.get("error")}', 'unresolved (sunucu çöküşü/zaman aşımı)')
        rec['status'] = 'FAIL'
        return rec, skrec, contrec
    if code == 400:
        msg = sheet.get('error', '')
        if verdict == 'other':
            add('geçerli girdi reddedildi', 'FAIL', msg, 'input validation (gereksiz ret)')
        elif verdict == 'discrete':
            consistent = ref['status'] == 'SINGULAR' or (ref['status'] == 'OK' and ref['cond'] > 1e10)
            add('ayrık teşhis ↔ scikit-fem', 'PASS' if consistent else 'FAIL',
                f'FEM Lab: {msg[:140]} · scikit-fem: {ref["status"]}, κ = {ref.get("cond", float("nan")):.2e}',
                None if consistent else 'numerical conditioning / diagnosis')
        add('ayrık sistem (bilgi)', 'PASS', f'scikit-fem: {ref["status"]}, κ = {ref.get("cond", float("nan")):.2e} '
            '(ayrık matrisin tekil olup olmaması sürekli hükmü kanıtlamaz)')
        add('tarayıcı çekirdeği aynı girdiyi', 'PASS' if js_kind == 'rejected' else 'WARN',
            js_val if js_kind != 'ok' else 'çözdü (tarayıcı sürekli düzeyi bilmez; ayrık sistem tekil değil)')
        statuses = [c['status'] for c in rec['checks']]
        rec['status'] = 'FAIL' if 'FAIL' in statuses else 'DIAGNOSED'
        return rec, skrec, contrec

    nodal = np.array(sheet['nodal'], dtype=float)
    text = sheet['text']
    items = {it['name']: it for it in sheet.get('verification', [])}
    cond = ref.get('cond') if ref['status'] == 'OK' else None
    add('kâğıt bütünlüğü', 'PASS' if ('nan' not in text.lower().replace('sınan', '') and 'None' not in text and '\\mathbf{a} =' in text)
        else 'FAIL', 'NaN/None yok, a vektörü yazılı', None)
    for name, it in items.items():
        if it['status'] == 'failed':
            add(f'uygulama doğrulaması: {name}', 'FAIL', it['detail'][:240], 'independent verification (uygulama Katman 4)')
        elif it['status'] == 'warning':
            add(f'uygulama doğrulaması: {name}', 'WARN', it['detail'][:240])
    ifem_item = next((it for n, it in items.items() if n.startswith('independent re-assembly')), None)
    M['FEM Lab ↔ ifem'] = APP_ITEM.get(ifem_item['status'], 'NOT_RUN') if ifem_item else 'NOT_RUN'
    fixed = {int(k): v for k, v in sheet.get('fixed', {}).items()}
    if fixed:
        bc_gap = max(abs(nodal[i] - v) for i, v in fixed.items())
        M['sınır koşulu'] = 'PASS' if bc_gap <= 1e-12 * max([1.0] + [abs(v) for v in fixed.values()]) else 'FAIL'
        add('temel koşul düğümde', M['sınır koşulu'], f'fark {bc_gap:.1e}', 'boundary-condition elimination')
    if js_kind == 'ok':
        tol = scale_tol(1e-9, cond)                 # tarayıcı: 96 noktalı bileşik Gauss; düzgün veride ~1e-15, 1e-9 geniş pay
        gap = rel_gap(js_val, nodal)
        M['tarayıcı ↔ sunucu'] = 'PASS' if gap <= tol else 'FAIL'
        add('tarayıcı ↔ SymPy düğüm', M['tarayıcı ↔ sunucu'], f'göreli fark {gap:.1e} (tolerans {tol:.0e}: 1e-9 ya da 100κε)',
            'element integration / assembly (iki çekirdekten biri)')
    else:
        M['tarayıcı ↔ sunucu'] = 'FAIL'
        add('tarayıcı ↔ SymPy düğüm', 'FAIL', f'tarayıcı {js_kind}: {js_val}', 'unresolved (tarayıcı çekirdeği)')
    if ref['status'] == 'OK':
        qtol = 1e-10 + 10 * ref['quadrature_estimate']
        if sheet.get('element_K'):
            M['scikit-fem eleman'], ev, worst = element_comparison(sheet['element_K'], sheet['element_F'], ref, qtol)
            add('scikit-fem ↔ FEM Lab eleman matrisleri ve yük vektörleri', M['scikit-fem eleman'], ev, 'element matrix construction')
            skrec['element_gap'] = worst
        if sheet.get('K'):
            gK = rel_gap(ref['K'], sheet['K'])
            M['scikit-fem global K'] = 'PASS' if gK <= qtol else 'FAIL'
            add('scikit-fem ↔ FEM Lab global K', M['scikit-fem global K'], f'göreli {gK:.1e} (tolerans {qtol:.0e}: 1e-10 + 10×kuadratür tahmini)',
                'global assembly')
            skrec['global_K_gap'] = gK
        else:
            M['scikit-fem global K'] = 'NOT_RUN'
            add('scikit-fem ↔ FEM Lab global K', 'NOT_RUN', 'K kâğıt yükünde yok (> 65 bilinmeyen)')
        gA = rel_gap(ref['a'], nodal)
        atol = max(qtol, scale_tol(1e-10, cond))
        M['scikit-fem düğüm'] = 'PASS' if gA <= atol else 'FAIL'
        add('scikit-fem ↔ FEM Lab düğüm katsayıları', M['scikit-fem düğüm'], f'göreli {gA:.1e} (tolerans {atol:.0e})',
            'unresolved (sistem ya da sınır koşulu)')
        skrec['nodal_gap'] = gA
        if ref['solver_residual'] > 1e-8 * max(1.0, float(np.max(np.abs(ref['a'])))) * max(1.0, (cond or 1) * EPS * 1e3):
            add('scikit-fem çözücü kalıntısı', 'WARN', f'{ref["solver_residual"]:.1e}')
        if cond and cond > 1e10:
            add('koşul sayısı', 'WARN', f'κ(K_red) = {cond:.2e} (scikit-fem, 2-norm)')
    elif ref['status'] == 'SINGULAR':
        for k in ('scikit-fem eleman', 'scikit-fem global K', 'scikit-fem düğüm'):
            M[k] = 'FAIL'
        add('scikit-fem', 'FAIL', f'referans sistem tekil (rank {ref["rank"]}/{ref["size"]}) ama FEM Lab çözdü', 'numerical conditioning / diagnosis')
    else:
        for k in ('scikit-fem eleman', 'scikit-fem global K', 'scikit-fem düğüm'):
            M[k] = 'UNSUPPORTED'
        add('scikit-fem', 'UNSUPPORTED', ref.get('reason', ref['status']))
    if sheet.get('exact'):
        valid, evidence = check_exact(inp, sheet['exact'])
        M['analitik çözüm'] = 'PASS' if valid else 'FAIL' if valid is False else 'INCONCLUSIVE'
        add('analitik çözüm (bağımsız sınama)', M['analitik çözüm'], evidence, None if valid else 'independent verification (analitik çözüm)')
        if valid:
            l2, h1 = error_norms(inp, nodal, sheet['exact'])
            s_l2 = sheet_number(text, r'\\\|u - u_h\\\|_\{L\^2\} = .*?\\approx ([-\d.]+(?:\\times 10\^\{-?\d+\})?)')
            s_h1 = sheet_number(text, r'\|u - u_h\|_\{H\^1\} = .*?\\approx ([-\d.]+(?:\\times 10\^\{-?\d+\})?)')
            norm_status = []
            for name, mine, theirs in (('L²', l2, s_l2), ('H¹', h1, s_h1)):
                if theirs is None:
                    norm_status.append('NOT_RUN')
                    add(f'kâğıttaki {name} hatası', 'NOT_RUN', 'kâğıtta bulunamadı')
                else:
                    ok = abs(mine - theirs) <= 1e-5 * max(mine, theirs) + 1e-13      # kâğıt 6 anlamlı basamak yazar
                    norm_status.append('PASS' if ok else 'FAIL')
                    add(f'kâğıttaki {name} hatası ↔ bağımsız', 'PASS' if ok else 'FAIL', f'kâğıt {theirs:.6g}, bağımsız {mine:.6g}',
                        None if ok else 'error calculation')
            M['L²/H¹ yeniden hesap'] = 'FAIL' if 'FAIL' in norm_status else 'PASS' if norm_status == ['PASS', 'PASS'] else 'NOT_RUN'
    else:
        M['analitik çözüm'] = 'INCONCLUSIVE'
        add('analitik çözüm', 'INCONCLUSIVE', sheet.get('exact_reason') or 'yok')
    rec['nodal'] = nodal.tolist()
    statuses = [c['status'] for c in rec['checks']]
    rec['status'] = 'FAIL' if 'FAIL' in statuses else 'WARN' if 'WARN' in statuses else 'PASS'
    if rec['status'] == 'WARN' and meta['validity'] == 'challenging':
        rec['status'] = 'LIMITATION'
    return rec, skrec, contrec


# Taban karşılaştırması: aynı girdinin değişiklik öncesi sunucu çıktısıyla (baseline_v1/server_outputs.json).
def baseline_comparison(results, base):
    rows = []
    for r in results:
        b = base.get(r['id'])
        if not b:
            continue
        before, after = b['code'], r['server']['code']
        num = None
        if before == 200 and after == 200 and r.get('nodal') is not None:
            num = rel_gap(r['nodal'], b['nodal'])
        rows.append({'id': r['id'], 'before': before, 'after': after, 'nodal_change': num,
                     'error_before': (b.get('error') or '')[:90], 'error_after': (r['server'].get('error') or '')[:90]})
    return rows


# Rapor: kapsam, yürütme, çapraz doğrulama matrisi, başarısızlıklar, scikit-fem (eleman/global/düğüm), sürekli düzey, taban farkı.
def make_report(corpus, results, skres, contres, base_rows, version):
    from collections import Counter
    L = []
    w = L.append
    w('# FEM Lab — 80 vakalık girdi korpusu: yürütme raporu\n')
    w(f'scikit-fem {version} · korpus: `verification/corpus/cases.json` (yalnız girdi, değişmedi) · sonuçlar: `results.json`, '
      '`skfem_results.json`, `continuous_results.json` · taban: `baseline_v1/`\n')
    w('## Kapsam\n')
    for title, key in (('Kategori', lambda c: c['meta']['category']), ('Geçerlilik', lambda c: c['meta']['validity']),
                       ('Mertebe', lambda c: '2' if (c['input'].get('alpha') or '').strip() not in ('', '0') else '1'),
                       ('Derece p', lambda c: str(c['input'].get('degree'))), ('n', lambda c: str(c['input'].get('n')))):
        cnt = Counter(key(c) for c in corpus)
        w(f'- **{title}:** ' + ', '.join(f'{k}: {v}' for k, v in sorted(cnt.items(), key=lambda kv: str(kv[0]))))
    cnt = Counter(c['input'][s]['type'] for c in corpus for s in ('left', 'right'))
    w('- **Koşul türleri (uç başına):** ' + ', '.join(f'{k}: {v}' for k, v in sorted(cnt.items())))
    w('\n## Yürütme özeti\n')
    st = Counter(r['status'] for r in results)
    w('| Durum | Vaka |\n|---|---:|')
    for k in ('PASS', 'WARN', 'LIMITATION', 'DIAGNOSED', 'REJECTED_OK', 'FAIL', 'RUNNER_ERROR'):
        w(f'| {k} | {st.get(k, 0)} |')
    w('\n## Çapraz doğrulama matrisi (geçerli + zorlayıcı vakalar)\n')
    w('Hücre değerleri: PASS · FAIL · WARN · UNSUPPORTED · INCONCLUSIVE · NOT_RUN. INCONCLUSIVE / UNSUPPORTED / NOT_RUN hiçbir '
      'zaman doğrulanmış sayılmaz.\n')
    short = {'PASS': '✓', 'FAIL': '✗', 'WARN': '!', 'UNSUPPORTED': 'U', 'INCONCLUSIVE': '?', 'NOT_RUN': '–'}
    w('| ID | Durum | ' + ' | '.join(MATRIX) + ' |\n|---|---|' + '---|' * len(MATRIX))
    tally = {k: Counter() for k in MATRIX}
    for r in results:
        if r['validity'] == 'invalid':
            continue
        for k in MATRIX:
            tally[k][r['matrix'][k]] += 1
        w(f"| {r['id']} | {r['status']} | " + ' | '.join(short[r['matrix'][k]] for k in MATRIX) + ' |')
    w('\n| Kontrol | ' + ' | '.join(short) + ' |\n|---|' + '---:|' * len(short))
    for k in MATRIX:
        w(f'| {k} | ' + ' | '.join(str(tally[k].get(s, 0)) for s in short) + ' |')
    w('\n## Başarısızlıklar (kanıt ve kaydedilen katman)\n')
    fails = [(r, c) for r in results for c in r['checks'] if c['status'] == 'FAIL']
    w('Yok.' if not fails else '\n'.join(f"- **{r['id']}** ({r['name']}) — {c['check']}: {c['evidence']} · katman: {c.get('layer', 'belirlenemedi')}"
                                          for r, c in fails))
    w('\n## Uyarılar\n')
    warns = [(r, c) for r in results for c in r['checks'] if c['status'] == 'WARN']
    w('Yok.' if not warns else '\n'.join(f"- {r['id']} — {c['check']}: {c['evidence']}" for r, c in warns))
    w('\n## scikit-fem bağımsız doğrulaması\n')
    sk = Counter(x['status'] for x in skres)
    w(f"- Çalıştırılan: {len(skres)} · çözüldü: {sk.get('OK', 0)} · tekil: {sk.get('SINGULAR', 0)} · desteklenmeyen: {sk.get('UNSUPPORTED', 0)}")
    el = [x for x in skres if x.get('element_gap')]
    if el:
        worst = max(el, key=lambda x: x['element_gap']['rel'])
        w(f"- Eleman düzeyi karşılaştırılan vaka: {len(el)}; en büyük göreli fark {worst['element_gap']['rel']:.2e} "
          f"(mutlak {worst['element_gap']['abs']:.2e}, {worst['id']} {worst['element_gap']['where']})")
    gk = [x['global_K_gap'] for x in skres if x.get('global_K_gap') is not None]
    na = [x['nodal_gap'] for x in skres if x.get('nodal_gap') is not None]
    if gk:
        w(f'- Global K karşılaştırılan: {len(gk)}; en büyük göreli fark {max(gk):.2e}')
    if na:
        w(f'- Düğüm katsayısı karşılaştırılan: {len(na)}; en büyük göreli fark {max(na):.2e}')
    w(f"- En büyük kuadratür tahmini (intorder 2p+40 ↔ 2p+8): {max([x['quadrature_estimate'] or 0 for x in skres if x.get('quadrature_estimate') is not None] + [0]):.2e}")
    w('- Taban dönüşümü: P1/P2 düğümsel (yalnız serbestlik derecesi sırası farklı); ElementLinePp(p ≥ 3) iç serbestlik derecelerinde '
      'hiyerarşik. Φ[J, m] = ψ_m(x_J), K_Lagrange = Φ⁻ᵀ K Φ⁻¹, F_Lagrange = Φ⁻ᵀ F (eleman ve global düzeyde).')
    w('\n## Sürekli problem doğrulaması\n')
    w('| ID | Sınıf | Düzenlilik (R) | Varlık | Teklik | Uyumluluk |\n|---|---|---|---|---|---|')
    for c in contres:
        reg = c.get('regularity', {}).get('established')
        w(f"| {c['id']} | {c.get('class', '?')} | {'evet' if reg else 'hayır' if reg is False else '?'} | {c['existence']['status']} | "
          f"{c['uniqueness']['status']} | {c['compatibility']['status']} |")
    cs = Counter(c['existence']['status'] for c in contres)
    w('\nVarlık durumları: ' + ', '.join(f'{k}: {v}' for k, v in sorted(cs.items())))
    w('\n### C64–C70 ayrıntıları\n')
    for c in contres:
        if c['id'] in ('C64', 'C65', 'C66', 'C67', 'C68', 'C69', 'C70'):
            w(f"**{c['id']}** ({c.get('class', '?')}; düzenlilik: {c.get('regularity', {}).get('evidence', '?')})")
            for part in ('existence', 'uniqueness', 'compatibility'):
                x = c[part]
                w(f"- {part}: **{x['status']}** — {x.get('method', '')} {x.get('evidence', '')}" + (f" ({x['cross_check']})" if x.get('cross_check') else ''))
            w('')
    w('## Taban (değişiklik öncesi) ile karşılaştırma\n')
    changed = [b for b in base_rows if b['before'] != b['after'] or (b['nodal_change'] or 0) > 1e-13]
    num = [b['nodal_change'] for b in base_rows if b['nodal_change'] is not None]
    w(f"- Her iki çalıştırmada kâğıt üretilen vaka: {len(num)}; düğüm değerlerinde en büyük göreli değişim: {max(num + [0]):.2e}")
    w('- HTTP durumu ya da düğüm değeri değişen vakalar:' + ('' if changed else ' yok'))
    for b in changed:
        w(f"  - {b['id']}: {b['before']} → {b['after']}" + (f", düğüm değişimi {b['nodal_change']:.1e}" if b['nodal_change'] else '')
          + (f" · önce: {b['error_before']}" if b['error_before'] else ''))
    return '\n'.join(L) + '\n'


def main(selected, jobs=4):
    from concurrent.futures import ThreadPoolExecutor
    corpus = json.loads(CORPUS.read_text(encoding='utf-8'))['cases']
    if selected:
        corpus = [c for c in corpus if c['id'] in selected]
    OUT.mkdir(exist_ok=True)

    def safe(case):
        try:
            return run_case(case)
        except Exception:                           # çalıştırıcının kendi hatası uygulama hatası sayılmaz; ayrı durum
            import traceback
            return ({'id': case['id'], 'name': case['meta']['name'], 'category': case['meta']['category'],
                     'validity': case['meta']['validity'], 'status': 'RUNNER_ERROR', 'server': {'code': None, 'seconds': None},
                     'matrix': {k: 'NOT_RUN' for k in MATRIX},
                     'checks': [{'check': 'çalıştırıcı', 'status': 'FAIL', 'evidence': traceback.format_exc()[-400:],
                                 'layer': 'test çalıştırıcısı (uygulama değil)'}]}, None, None)
    with ThreadPoolExecutor(jobs) as ex:
        outs = list(ex.map(safe, corpus))
    results = [o[0] for o in outs]
    skres = [o[1] for o in outs if o[1]]
    contres = [o[2] for o in outs if o[2]]
    for r in results:
        print(f"{r['id']} {r['status']:<12} {r['category']:<28} {r['name'][:60]}", flush=True)
    import skfem
    base_file = OUT / 'baseline_v1' / 'server_outputs.json'
    base_rows = baseline_comparison(results, json.loads(base_file.read_text(encoding='utf-8'))) if base_file.exists() else []
    dump = lambda name, obj: (OUT / name).write_text(json.dumps(obj, ensure_ascii=False, indent=1, default=str), encoding='utf-8')
    dump('results.json', {'skfem_version': skfem.__version__, 'results': results})
    dump('skfem_results.json', {'skfem_version': skfem.__version__, 'results': skres})
    dump('continuous_results.json', {'results': contres})
    dump('baseline_comparison.json', {'rows': base_rows})
    (OUT / 'report.md').write_text(make_report(json.loads(CORPUS.read_text(encoding='utf-8'))['cases'] if not selected else corpus,
                                               results, skres, contres, base_rows, skfem.__version__), encoding='utf-8')
    return results


if __name__ == '__main__':
    res = main(sys.argv[1:])
    bad = [r for r in res if r['status'] in ('FAIL', 'RUNNER_ERROR')]
    print(f'\n{len(res)} vaka · FAIL {len(bad)}')
    sys.exit(1 if bad else 0)
