#!/usr/bin/env python3
"""run_verification.py — bağımsız doğrulama koşucusu.

    python3 verification/run_verification.py --quick        (n = 4, 8; ~1 dk)
    python3 verification/run_verification.py --full         (n = 4…32; birkaç dk)
    python3 verification/run_verification.py --baseline     (temel düğüm değerlerini kaydeder)
    python3 verification/run_verification.py --compare      (temel kayıtla toleranslı karşılaştırma)

Neden: uygulamanın sonuçları (tarayıcı çekirdeği, `node tests.js --solve`) bağımsız FEM'e (ifem.py) ve analitik çözüme
karşı ölçülür. Beklenen değerler uygulamanın kendi çıktısından alınmaz. Çıkış kodu 0 = hepsi geçti.
"""
import json, os, subprocess, sys
import mpmath as mp

HERE = os.path.dirname(os.path.abspath(__file__))
APP = os.path.dirname(HERE)
sys.path.insert(0, HERE)
import ifem, benchmarks as B   # noqa: E402

AGREE_REL = mp.mpf('1e-9')        # uygulama ↔ bağımsız katsayı farkı / max|c| (iki 64-bit/30-basamak çözüm arası; sezgisel)
FLOOR = mp.mpf('1e-11')           # bu L² altı "tabanda" sayılır, hız hesaplanmaz (uygulamadaki taban ile aynı mertebe)
BASE_REL = 1e-12                  # temel kayıtla regresyon toleransı (göreli, max|c| ölçekli)


# Uygulamanın tarayıcı çekirdeğiyle çözüm (aynı yol tests.py'de de kullanılır).
def app_solve(spec, n, p):
    run = subprocess.run(['node', os.path.join(APP, 'tests.js'), '--solve', json.dumps({**spec, 'n': n, 'degree': p})],
                         capture_output=True, text=True)
    if run.returncode:
        return None, (run.stderr.strip().splitlines() or ['?'])[-1]
    return json.loads(run.stdout), None


# Tek benchmark: her p için n listesinde bağımsız ↔ uygulama farkı, hata ölçüleri ve son geçerli çiftten L² hızı.
def run_case(case, ns):
    code, cls, spec, ind, u, du, expect, note = case
    lines, ok = [], True
    for p in (1, 2, 3, 4):
        prev, rate, worst = None, None, mp.mpf(0)
        for n in ns:
            c_app, err = app_solve(spec, n, p)
            if c_app is None:
                lines.append(f'  p={p} n={n}: UYGULAMA HATASI {err}'); ok = False; break
            s = ifem.solve(ind, n, p)
            scale = max(1, max(abs(v) for v in s['c']))
            worst = max(worst, max(abs(mp.mpf(a) - b) for a, b in zip(c_app, s['c'])) / scale)
            e = ifem.errors(s, u, du, c_override=c_app)
            if prev is not None and prev > FLOOR and e['l2'] > FLOOR:
                rate = mp.log(prev / e['l2'], 2)
            prev = e['l2']
        want = expect[p - 1]
        if want not in (None, 'exact') and ns[-2] < B.ASYMPTOTIC_FROM.get(code, 0):
            want = None                                  # bu ağlar ön-asimptotik: hız iddiası yok (rapor edilir)
        agree = worst <= AGREE_REL
        if want == 'exact':
            rate_ok = prev is not None and prev <= FLOOR * 100
        elif want is None:
            rate_ok = True
        else:
            rate_ok = rate is not None and rate >= want
        ok &= agree and rate_ok
        lines.append(f'  p={p}: L2(son)={mp.nstr(prev, 3)} hız={mp.nstr(rate, 3) if rate is not None else "—"} '
                     f'beklenen≥{want} | app−bağımsız göreli {mp.nstr(worst, 2)} {"✓" if agree and rate_ok else "✗"}')
    return ok, lines


# Temel kayıt: tüm benchmark'larda tarayıcı çekirdeğinin düğüm değerleri (regresyon karşılaştırması için).
def baseline(path, ns):
    data = {}
    for case in B.CASES:
        for p in (1, 2, 3, 4):
            for n in ns:
                c, err = app_solve(case[2], n, p)
                data[f'{case[0]}|p{p}|n{n}'] = c if c is not None else {'error': err}
    json.dump(data, open(path, 'w'), indent=0)
    return data


def compare(path, ns):
    old = json.load(open(path))
    worst, bad = 0.0, []
    for key, c_old in old.items():
        code, p, n = key.split('|')
        case = next(c for c in B.CASES if c[0] == code)
        c_new, err = app_solve(case[2], int(n[1:]), int(p[1:]))
        if isinstance(c_old, dict) or c_new is None:
            if not (isinstance(c_old, dict) and c_new is None):
                bad.append(f'{key}: davranış değişti ({c_old if isinstance(c_old, dict) else "çözüm"} → {err or "çözüm"})')
            continue
        scale = max(1.0, max(abs(v) for v in c_old))
        d = max(abs(a - b) for a, b in zip(c_old, c_new)) / scale
        worst = max(worst, d)
        if d > BASE_REL:
            bad.append(f'{key}: göreli fark {d:.2e}')
    return worst, bad


if __name__ == '__main__':
    mode = sys.argv[1] if len(sys.argv) > 1 else '--quick'
    base = os.path.join(HERE, 'baseline', 'nodal_js.json')
    if mode == '--baseline':
        d = baseline(base, [4, 8, 16, 32]); print(f'{len(d)} kayıt → {base}'); sys.exit(0)
    if mode == '--compare':
        worst, bad = compare(base, None)
        print(f'temel kayıtla en büyük göreli fark: {worst:.2e} (tolerans {BASE_REL:.0e}); sorunlu: {len(bad)}')
        for b in bad[:20]: print('  ', b)
        sys.exit(1 if bad else 0)
    ns = [4, 8] if mode == '--quick' else [4, 8, 16, 32]
    failures = 0
    for case in B.CASES:
        ok, lines = run_case(case, ns)
        print(f'{"OK  " if ok else "FAIL"} {case[0]} [{case[1]}] {case[7]}')
        for l in lines: print(l)
        failures += not ok
    print(f'\n{len(B.CASES) - failures}/{len(B.CASES)} benchmark geçti ({mode})')
    sys.exit(1 if failures else 0)
