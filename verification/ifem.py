"""ifem.py — BAĞIMSIZ doğrulama FEM'i (uygulamanın hiçbir dosyasını içe aktarmaz).

Neden var: uygulamanın iki çekirdeği (index.html JS, fem_sheet.py SymPy) birbirini doğrular ama aynı kişinin aynı
kararlarını taşır. Bu dosya aynı formülasyonu FARKLI bir gerçekleştirmeyle çözer:
  - şekil fonksiyonları ürün formülü yerine Vandermonde tersinden (N_k = Σ_m (V⁻¹)_{mk} t^m),
  - kuadratür 24 noktalı Gauss–Legendre (mpmath, 30 basamak),
  - çözücü mpmath LU,
  - katsayılar ayrıştırıcıdan geçmez, doğrudan Python fonksiyonudur.
Formülasyon (uygulamanın beyan ettiği): 2. mertebe ∫(v'αu' − vβu' − vγu) = [vαu'] − ∫vf;
1. mertebe ∫v(βu' + γu) = ∫vf; u verilen düğümün satırı atılır; IVP'de akısı bilinmeyen uçtaki satır atılır.
Bu, kod bağımsızlığıdır; matematiksel kanıt değildir.
"""
import mpmath as mp

mp.mp.dps = 30
_GL = mp.calculus.quadrature.GaussLegendre(mp.mp)
NODES = [(mp.mpf(t), mp.mpf(w)) for t, w in _GL.calc_nodes(4, mp.mp.prec)]   # 24 nokta, [-1, 1]


# Bant LU. Neden: mpmath.lu_solve / A**-1 yoğun matriste O(n³) saf Python; 129–257 bilinmeyende doğrulama tek başına
# 10–15 s, CI makinesinde paralel yük altında sunucunun 60 s sınırını aştı (korpus C48/C50, HTTP 504). Burada mpmath.LU_decomp,
# L_solve, U_solve ve inverse'ün AYNI aritmetiği uygulanır (prec + 10, satır toplamıyla ölçekli kısmi pivot, aynı tekillik
# eşiği ve aynı işlem sırası); yalnız tam sıfır olan bant dışı elemanlarla yapılan işlemler atlanır (bunlar değeri
# değiştirmez). Uygulamanın fem_sheet.band_* fonksiyonları ÇAĞRILMAZ: doğrulama katmanı motordan bağımsız kalır.
def _band_lu(A):
    n = len(A)
    kl = max([i - j for i, row in enumerate(A) for j, v in enumerate(row) if v != 0] + [0])
    ku = max([j - i for i, row in enumerate(A) for j, v in enumerate(row) if v != 0] + [0])
    reach = kl + ku                                   # kısmi pivot sonrası U'nun üst bant genişliği en çok kl + ku
    A = [list(row) for row in A]
    tol = mp.absmin(mp.mnorm(mp.matrix(A), 1) * mp.eps)
    # j = 0 adımında mpmath bütün satırların |·| toplamını sınar; bandın altında kalan satırlara eleme dokunmadığı için
    # toplamları sonraki adımlarda değişmez — bir kez sınamak mpmath'in bütün adımlardaki sınamasına denktir.
    if any(mp.absmin(mp.fsum([mp.absmin(v) for v in row])) <= tol for row in A):
        raise ZeroDivisionError('matrix is numerically singular')
    p = [None] * (n - 1)
    for j in range(n - 1):
        biggest, last = 0, min(n - 1, j + kl)
        for k in range(j, last + 1):
            s = mp.fsum([mp.absmin(A[k][l]) for l in range(j, min(n, j + reach + kl + 1))])
            if mp.absmin(s) <= tol:
                raise ZeroDivisionError('matrix is numerically singular')
            current = 1 / s * mp.absmin(A[k][j])
            if current > biggest:
                biggest, p[j] = current, k
        A[j], A[p[j]] = A[p[j]], A[j]
        if mp.absmin(A[j][j]) <= tol:
            raise ZeroDivisionError('matrix is numerically singular')
        for i in range(j + 1, last + 1):
            A[i][j] /= A[j][j]
            for k in range(j + 1, min(n, j + reach + 1)):
                A[i][k] -= A[i][j] * A[j][k]
    if mp.absmin(A[n - 1][n - 1]) <= tol:
        raise ZeroDivisionError('matrix is numerically singular')
    # Satır değişimleri eski L çarpanlarını da taşıdığından L'nin bandı önceden bilinmez: sıfır olmayanlar listelenir.
    lower = [[(j, A[i][j]) for j in range(i) if A[i][j] != 0] for i in range(n)]
    upper = [[(j, A[i][j]) for j in range(i + 1, n) if A[i][j] != 0] for i in range(n)]
    return {'n': n, 'p': p, 'lower': lower, 'upper': upper, 'diag': [A[i][i] for i in range(n)]}


# mpmath L_solve + U_solve (aynı işlem sırası); çağıran prec + 10'da tutar.
def _band_substitute(lu, b):
    b = list(b)
    for k, pk in enumerate(lu['p']):
        b[k], b[pk] = b[pk], b[k]
    for i in range(1, lu['n']):
        for j, v in lu['lower'][i]:
            b[i] -= v * b[j]
    x = list(b)
    for i in range(lu['n'] - 1, -1, -1):
        for j, v in lu['upper'][i]:
            x[i] -= v * x[j]
        x[i] /= lu['diag'][i]
    return x


# mpmath.lu_solve(A, b) ile aynı sonuç (kare A). A liste listesi, b liste; dönüş mp.matrix sütunu.
def band_lu_solve(A, b):
    with mp.workprec(mp.mp.prec + 10):
        lu = _band_lu(A)
        x = _band_substitute(lu, [mp.mpf(v) for v in b])
    return mp.matrix(x)


# mpmath.inverse(A) (= A**-1) ile aynı sonuç: aynı LU, birim vektörlerle bant çözümü. Dönüş mp.matrix.
def band_inverse(A):
    n = len(A)
    with mp.workprec(mp.mp.prec + 10):
        lu = _band_lu(A)
        cols = [_band_substitute(lu, [mp.mpf(1) if k == i else mp.mpf(0) for k in range(n)]) for i in range(n)]
        return mp.matrix([[cols[j][i] for j in range(n)] for i in range(n)])


# Şekil fonksiyonları polinom katsayılarından: V c = â sisteminin tersi. Uygulamanın ürün formülünden bağımsız yol.
def basis(p):
    ts = [mp.mpf(j) / p for j in range(p + 1)]
    V = mp.matrix([[t ** m for m in range(p + 1)] for t in ts])
    C = V ** -1
    # Değerler elemandan bağımsızdır (yalnız t ve hassasiyete bağlı): her eleman ve her hata noktasında aynı t için yeniden
    # hesaplamak yerine bir kez hesaplanıp saklanır. Aritmetik ve sonuç aynıdır; kopya döner (çağıran listeyi paylaşmaz).
    # Neden: n = 64, p = 4'te 2.9 milyon t**m çağrısı doğrulama süresinin çoğuydu (korpus C50, CI süre sınırı).
    memo_N, memo_dN = {}, {}

    def N(t):
        key = (t, mp.mp.prec)
        if key not in memo_N:
            memo_N[key] = [sum(C[m, k] * t ** m for m in range(p + 1)) for k in range(p + 1)]
        return list(memo_N[key])

    def dN(t):
        key = (t, mp.mp.prec)
        if key not in memo_dN:
            memo_dN[key] = [sum(C[m, k] * m * t ** (m - 1) for m in range(1, p + 1)) for k in range(p + 1)]
        return list(memo_dN[key])
    return ts, N, dN, C


# Montaj ve sınır koşulu uygulaması; uygulamanın beyan ettiği satır seçimi kuralıyla.
def solve(prob, n, p):
    a_, b_, g_, f_ = prob['alpha'], prob['beta'], prob['gamma'], prob['f']
    x0, x1 = mp.mpf(prob['x0']), mp.mpf(prob['x1'])
    h, dofs, order = (x1 - x0) / n, n * p + 1, prob['order']
    _, N, dN, _ = basis(p)
    K, F = mp.zeros(dofs, dofs), mp.zeros(dofs, 1)
    for e in range(n):
        xe = x0 + e * h
        for t, w in NODES:
            tt = (t + 1) / 2
            x, wx = xe + h * tt, w * h / 2
            Nv, dNv = N(tt), [d / h for d in dN(tt)]
            A, B, G, Fx = a_(x), b_(x), g_(x), f_(x)
            for i in range(p + 1):
                gi = e * p + i
                F[gi] += wx * (Nv[i] * Fx if order == 1 else -Nv[i] * Fx)
                for j in range(p + 1):
                    K[gi, e * p + j] += wx * ((Nv[i] * (B * dNv[j] + G * Nv[j])) if order == 1 else
                                              (dNv[i] * A * dNv[j] - Nv[i] * B * dNv[j] - Nv[i] * G * Nv[j]))
    L, R = prob['left'], prob['right']
    if order == 2:
        for side, idx, xb, sgn in (('left', 0, x0, -1), ('right', dofs - 1, x1, 1)):
            c = prob[side]
            if c[0] == 'neumann':
                F[idx] += sgn * a_(xb) * c[1]
            if c[0] == 'flux':                     # akı koşulu α u' = q: Q = q (türev koşuluna çevrilmez)
                F[idx] += sgn * c[1]
            if c[0] == 'both':
                F[idx] += sgn * a_(xb) * c[2]
    fixed = {i: mp.mpf(prob[s][1]) for s, i in (('left', 0), ('right', dofs - 1)) if prob[s][0] in ('dirichlet', 'both')}
    free = [i for i in range(dofs) if i not in fixed]
    rows = list(range(dofs - 1)) if L[0] == 'both' else list(range(1, dofs)) if R[0] == 'both' else free
    if not free:                       # tüm düğümler temel koşulla belirli: çözülecek sistem yok (boş matris kurulmaz)
        sol = []
    else:
        A = [[K[i, j] for j in free] for i in rows]
        bb = [F[i] - sum(K[i, j] * v for j, v in fixed.items()) for i in rows]
        sol = band_lu_solve(A, bb)              # mp.lu_solve ile aynı aritmetik, bant içinde (bkz. _band_lu)
    c = [None] * dofs
    for i, v in fixed.items():
        c[i] = v
    for k, i in enumerate(free):
        c[i] = sol[k]
    return {'c': c, 'K': K, 'F': F, 'h': h, 'n': n, 'p': p, 'x0': x0, 'free': free, 'rows': rows, 'fixed': fixed}


# Hata ölçüleri: L², H¹-yarı, düğüm ve E_max. E_max: eleman başına 200 eşit örnek + en büyüğün çevresinde
# türevin sıfırı (findroot). Uygulamanın yönteminden (32 nokta + altın oran) farklı bir yol.
def errors(sol, u, du, c_override=None):
    s = dict(sol)
    if c_override is not None:
        s['c'] = [mp.mpf(v) for v in c_override]
    n, p, h, x0 = s['n'], s['p'], s['h'], s['x0']
    _, N, dN, _ = basis(p)
    l2 = h1 = emax = nodal = mp.mpf(0)
    for e in range(n):
        xe = x0 + e * h
        uh = lambda t, e=e: sum(s['c'][e * p + k] * v for k, v in enumerate(N(t)))
        duh = lambda t, e=e: sum(s['c'][e * p + k] * v / h for k, v in enumerate(dN(t)))
        for t, w in NODES:
            tt = (t + 1) / 2
            x, wx = xe + h * tt, w * h / 2
            l2 += wx * (u(x) - uh(tt)) ** 2
            h1 += wx * (du(x) - duh(tt)) ** 2
        ts = [mp.mpf(k) / 200 for k in range(201)]
        vals = [abs(u(xe + h * t) - uh(t)) for t in ts]
        k = max(range(201), key=lambda i: vals[i])
        lo, hi = ts[max(0, k - 1)], ts[min(200, k + 1)]
        g = lambda t: -abs(u(xe + h * t) - uh(t))
        try:
            tm = mp.findroot(lambda t: mp.diff(g, t), ts[k])
            tm = tm if lo <= tm <= hi else ts[k]
        except Exception:
            tm = ts[k]
        emax = max(emax, vals[k], -g(tm))
    for i, cv in enumerate(s['c']):
        nodal = max(nodal, abs(u(x0 + i * h / p) - cv))
    return {'l2': mp.sqrt(l2), 'h1': mp.sqrt(h1), 'max': emax, 'nodal': nodal}
