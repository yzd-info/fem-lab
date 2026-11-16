"""skfem_reference.py — scikit-fem ile BAĞIMSIZ sonlu eleman referansı (FEM Lab'ın hiçbir dosyasını içe aktarmaz).

Neden var: FEM Lab'ın iki çekirdeği ve verification/ifem.py aynı geliştirme sürecinden çıktı. scikit-fem, ayrı bir ekibin
yazdığı yerleşik bir FEM kütüphanesidir; zayıf formu kendi derleyicisiyle (BilinearForm / LinearForm, kendi kuadratürü,
kendi eleman tabanı) montajlar. FEM Lab'ın montajlı matrisleri BURAYA VERİLMEZ: girdi metni okunur, zayıf form baştan kurulur.

Kullanılan formülasyon (FEM Lab'ın belgelediği problem tanımı — docs/MATEMATIK.md; ortak-mod hatası riski aşağıda):
  2. mertebe (α ≢ 0), deneme u_h, ağırlık v:
      a(u, v) = ∫ (α u' v' − β u' v − γ u v) dx,   ℓ(v) = −∫ f v dx,   F_b = [−Q₀, 0, …, 0, Q₁]ᵀ,   Q = α u'
      Q: türev koşulu u'(x_b) = g → Q = α(x_b)·g · akı koşulu α u'(x_b) = q → Q = q · u ve u' → Q = α(x_b)·s ·
         verilmedi (sınır değer) → Q = 0 · u verilen uç ve başlangıç değer probleminin serbest ucu → Q bilinmez.
  1. mertebe (α ≡ 0):  a(u, v) = ∫ (β u' + γ u) v dx,  ℓ(v) = ∫ f v dx  (kısmi integrasyon yok).
  Satır/sütun seçimi: u verilen düğümün sütunu bilinen değerle sağa taşınır. Kullanılan satırlar: sınır değer probleminde
  serbest düğümler; başlangıç değer probleminde akısı bilinmeyen uçtaki satır çıkarılır; 1. mertebede u verilen düğümün satırı.

Taban: ElementLineP1, ElementLineP2 (düğümsel Lagrange) ve ElementLinePp(p) (p ≥ 3; iç serbestlik dereceleri HİYERARŞİK,
düğümsel değil). Karşılaştırma için her şey FEM Lab'ın eşit aralıklı düğümsel (Lagrange) tabanına çevrilir:
  Φ[J, m] = ψ_m(x_J)  (scikit-fem tabanı ψ_m'nin FEM Lab düğümü x_J'deki değeri),  u_düğüm = Φ c,
  K_Lagrange = Φ⁻ᵀ K_skfem Φ⁻¹,  F_Lagrange = Φ⁻ᵀ F_skfem   (ψ_m = Σ_J Φ[J, m] N_J olduğu için).

Kuadratür: Gauss–Legendre, intorder = 2p + 40 (karşılaştırma) ve 2p + 8 (fark = kuadratür hatası tahmini). Katsayı ve yük
polinom değilse sonuçlar kuadratürle sınırlıdır; bu fark ayrıca raporlanır, montaj hatası sayılmaz.

Ortak-mod uyarısı: bu modül FEM Lab ile AYNI zayıf formu kullanır (problem tanımı gereği). İkisi aynı yanlış formülü
paylaşıyorsa uyuşabilirler; bu yüzden analitik çözümler, elle türetilmiş eleman matrisleri ve güçlü-form kalıntısı ayrıca sınanır.
"""
import numpy as np
import sympy as sp

import skfem
from skfem import MeshLine, Basis, BilinearForm, LinearForm, asm, ElementLineP1, ElementLineP2, ElementLinePp

from ref_parser import X, expr   # noqa: F401  ortak okuyucu (continuous_problem ile aynı)

Q_HIGH, Q_LOW = 40, 8          # intorder = 2p + Q_*
VERSION = skfem.__version__


def as_function(e):
    fn = sp.lambdify(X, e, 'numpy')
    return lambda x: np.broadcast_to(np.asarray(fn(x), dtype=float), np.shape(x)) * 1.0


def element_for(p):
    return ElementLineP1() if p == 1 else ElementLineP2() if p == 2 else ElementLinePp(p)


# Φ: scikit-fem serbestlik derecelerinden FEM Lab düğüm değerlerine geçiş (yerel taban elem.lbasis ile, referans t ∈ [0, 1]).
def change_of_basis(basis, elem, n, p):
    Phi = np.zeros((n * p + 1, basis.N))
    t = np.linspace(0, 1, p + 1)[None, :]
    local = np.array([elem.lbasis(t, i)[0] for i in range(p + 1)])      # (yerel dof i, yerel düğüm k)
    for e in range(n):
        for k in range(p + 1):
            for i in range(p + 1):
                Phi[e * p + k, basis.element_dofs[i, e]] = local[i, k]
    return Phi


def forms(order, A, Bt, G, Fs):
    if order == 2:
        a = BilinearForm(lambda u, v, w: A(w.x[0]) * u.grad[0] * v.grad[0] - Bt(w.x[0]) * u.grad[0] * v - G(w.x[0]) * u * v)
        l = LinearForm(lambda v, w: -Fs(w.x[0]) * v)
    else:
        a = BilinearForm(lambda u, v, w: Bt(w.x[0]) * u.grad[0] * v + G(w.x[0]) * u * v)
        l = LinearForm(lambda v, w: Fs(w.x[0]) * v)
    return a, l


def assemble(points, p, order, fns, intorder):
    mesh = MeshLine(np.asarray(points, dtype=float))
    elem = element_for(p)
    basis = Basis(mesh, elem, intorder=intorder)
    a, l = forms(order, *fns)
    K, F = asm(a, basis).toarray(), asm(l, basis)
    Phi = change_of_basis(basis, elem, len(points) - 1, p)
    Pinv = np.linalg.inv(Phi)
    # numpy 2.0 + macOS Accelerate matmul'da yanlış kayan nokta uyarıları verebiliyor; uyarı yerine sonuç açıkça denetlenir.
    with np.errstate(all='ignore'):
        K_lagr, F_lagr = Pinv.T @ K @ Pinv, Pinv.T @ F
    if not (np.all(np.isfinite(K_lagr)) and np.all(np.isfinite(F_lagr))):
        raise FloatingPointError('taban dönüşümünde sonlu olmayan değer')
    return K_lagr, F_lagr


# Ana giriş: girdi sözlüğü (FEM Lab /api/sheet ile aynı alanlar) → yapılandırılmış sonuç. Desteklenmeyen sınıf açıkça bildirilir.
def reference(inp):
    out = {'skfem_version': VERSION, 'status': 'NOT_RUN'}
    try:
        alpha, beta, gamma, f = (expr(inp.get(k)) for k in ('alpha', 'beta', 'gamma', 'f'))
        x0, x1 = float(expr(inp['x0'])), float(expr(inp['x1']))
        n, p = int(inp['n']), int(inp['degree'])
    except Exception as error:
        return {**out, 'status': 'UNSUPPORTED', 'reason': f'girdi okunamadı ({type(error).__name__})'}
    order = 2 if alpha != 0 else 1
    side = {s: (inp[s].get('type', 'none'), inp[s].get('value', ''), inp[s].get('slope', '')) for s in ('left', 'right')}
    ivp = any(side[s][0] == 'both' for s in side)
    A, Bt, G, Fs = (as_function(e) for e in (alpha, beta, gamma, f))
    nodes = np.linspace(x0, x1, n * p + 1)
    edges = np.linspace(x0, x1, n + 1)
    smooth = all(e.is_polynomial(X) for e in (alpha, beta, gamma, f))
    q_hi, q_lo = 2 * p + Q_HIGH, 2 * p + Q_LOW
    K, Fl = assemble(edges, p, order, (A, Bt, G, Fs), q_hi)
    K_lo, Fl_lo = assemble(edges, p, order, (A, Bt, G, Fs), q_lo)
    quad = max(np.max(np.abs(K - K_lo)) / max(1.0, np.max(np.abs(K))), np.max(np.abs(Fl - Fl_lo)) / max(1.0, np.max(np.abs(Fl))))
    elements = [assemble([edges[e], edges[e + 1]], p, order, (A, Bt, G, Fs), q_hi) for e in range(n)]

    N = n * p + 1
    Fb = np.zeros(N)
    Q = {}
    for s, xb, idx, sgn in (('left', x0, 0, -1), ('right', x1, N - 1, 1)):
        kind, value, slope = side[s]
        if order == 1:
            continue
        ab = float(A(np.array([xb]))[0])
        q = (ab * float(expr(value)) if kind == 'neumann' else float(expr(value)) if kind == 'flux'
             else ab * float(expr(slope)) if kind == 'both' else 0.0 if kind == 'none' and not ivp else None)
        Q[s] = q
        if q is not None:
            Fb[idx] += sgn * q
    fixed = {}
    for s, idx in (('left', 0), ('right', N - 1)):
        if side[s][0] in ('dirichlet', 'both'):
            fixed[idx] = float(expr(side[s][1]))
    free = [i for i in range(N) if i not in fixed]
    if side['left'][0] == 'both':
        rows = list(range(N - 1))
    elif side['right'][0] == 'both':
        rows = list(range(1, N))
    else:
        rows = free
    F = Fl + Fb
    Kred = K[np.ix_(rows, free)]
    Fred = F[rows] - sum(K[rows, j] * v for j, v in fixed.items()) if fixed else F[rows]
    cond = float(np.linalg.cond(Kred)) if Kred.size else 1.0
    a = np.zeros(N)
    for j, v in fixed.items():
        a[j] = v
    if Kred.size:
        try:
            sol = np.linalg.solve(Kred, Fred)
        except np.linalg.LinAlgError:
            return {**out, 'status': 'SINGULAR', 'cond': cond, 'K': K.tolist(), 'F_l': Fl.tolist(),
                    'rank': int(np.linalg.matrix_rank(Kred)), 'size': len(free)}
        for j, v in zip(free, sol):
            a[j] = v
        with np.errstate(all='ignore'):                # Accelerate'in yanlış uyarısı; sonuç aşağıda açıkça denetlenir
            residual = float(np.linalg.norm(Kred @ sol - Fred, np.inf))
        if not np.isfinite(residual):
            raise FloatingPointError('çözücü kalıntısı sonlu değil')
    else:
        residual = 0.0
    return {**out, 'status': 'OK', 'element': type(element_for(p)).__name__ + (f'({p})' if p >= 3 else ''),
            'degree': p, 'n': n, 'intorder': q_hi, 'intorder_check': q_lo, 'quadrature_estimate': float(quad),
            'polynomial_data': bool(smooth), 'formulation': ('order2-weak-ibp' if order == 2 else 'order1-direct-weighted-residual')
            + ('+ivp-petrov' if ivp else ''), 'nodes': nodes.tolist(), 'K': K.tolist(), 'F_l': Fl.tolist(), 'F_b': Fb.tolist(),
            'element_K': [k.tolist() for k, _ in elements], 'element_F': [g.tolist() for _, g in elements],
            'rows': rows, 'free': free, 'fixed': {str(k): v for k, v in fixed.items()}, 'Q': Q,
            'cond': cond, 'solver_residual': residual, 'a': a.tolist()}


# FEM Lab'ın düğüm katsayılarından (ya da bu modülün kendi çözümünden) u_h ve u_h' değerleri: eleman eleman Lagrange
# tabanı ile (FEM Lab'ın düğümsel tabanında yeniden kurulum; scikit-fem tabanı Φ ile buna eşdeğerdir).
def evaluate(a, x0, x1, n, p, points):
    a = np.asarray(a, dtype=float)
    h = (x1 - x0) / n
    t_nodes = np.linspace(0, 1, p + 1)
    vals, ders = [], []
    for x in points:
        e = min(int((x - x0) / h), n - 1)
        t = (x - x0) / h - e
        L = np.ones(p + 1); dL = np.zeros(p + 1)
        for k in range(p + 1):
            others = [m for m in range(p + 1) if m != k]
            L[k] = np.prod([(t - t_nodes[m]) / (t_nodes[k] - t_nodes[m]) for m in others])
            dL[k] = sum(np.prod([(t - t_nodes[m]) / (t_nodes[k] - t_nodes[m]) for m in others if m != j]) / (t_nodes[k] - t_nodes[j])
                        for j in others)
        coef = a[e * p:(e + 1) * p + 1]
        vals.append(float(coef @ L)); ders.append(float(coef @ dL / h))
    return np.array(vals), np.array(ders)
