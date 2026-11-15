"""benchmarks.py — analitik çözümü bilinen doğrulama problemleri (T1–T16) ve iyi tanımlılık sınır durumları (E1–E5).

Neden ayrı: beklenen değerler uygulamanın çıktısından DEĞİL, analitik çözümden ve problem matematiğinden gelir.
Her girdi: uygulama girdisi (metin), bağımsız girdi (Python fonksiyonları), u, u', beklenen hız (L², h-inceltme, asimptotik).
expected_rate: [p=1, p=2, p=3, p=4] için alt sınır (None = iddia yok; 'exact' = polinom çözüm, hata tabanda).
"""
import mpmath as mp

pi, sin, cos, exp, cosh, sinh, log = mp.pi, mp.sin, mp.cos, mp.exp, mp.cosh, mp.sinh, mp.log
z = lambda x: 0
D = lambda v: {'type': 'dirichlet', 'value': str(v)}
Nm = lambda v: {'type': 'neumann', 'value': str(v)}
BOTH = lambda u, s: {'type': 'both', 'value': str(u), 'slope': str(s)}
NONE = {'type': 'none', 'value': ''}


# Üretilmiş çözüm (T7) için f, sembolik türetmeyle elde edilmiş kapalı ifadedir: L u = f, u = sin πx + x².
def _manufactured():
    U = lambda x: sin(pi * x) + x ** 2
    dU = lambda x: pi * cos(pi * x) + 2 * x
    d2U = lambda x: -pi ** 2 * sin(pi * x) + 2
    f = lambda x: 2 * x * dU(x) + (1 + x ** 2) * d2U(x) + sin(x) * dU(x) + exp(x) * U(x)
    return U, dU, f


U7, dU7, f7 = _manufactured()
Bh = (1 - cos(10)) / sin(10)

CASES = [
    # (kod, sınıf, uygulama girdisi, bağımsız girdi, u, u', beklenen L² hızı alt sınırları, not)
    ('T1', '1. mertebe', dict(alpha='', beta='1', gamma='', f='2x', x0=0, x1=2, left=D(0), right=NONE),
     dict(alpha=z, beta=lambda x: 1, gamma=z, f=lambda x: 2 * x, x0=0, x1=2, left=('dirichlet', 0), right=('none',), order=1),
     lambda x: x ** 2, lambda x: 2 * x, [1.8, 'exact', 'exact', 'exact'], 'ödev kıyası'),
    ('T2', '1. mertebe', dict(alpha='', beta='1', gamma='', f='3x^2', x0=0, x1=1, left=D(0), right=NONE),
     dict(alpha=z, beta=lambda x: 1, gamma=z, f=lambda x: 3 * x ** 2, x0=0, x1=1, left=('dirichlet', 0), right=('none',), order=1),
     lambda x: x ** 3, lambda x: 3 * x ** 2, [1.8, 1.8, 'exact', 'exact'], '1. mertebede P2 hızı 2 (garanti yok)'),
    ('T3', '1. mertebe, reaksiyonlu', dict(alpha='', beta='1', gamma='2', f='0', x0=0, x1=1, left=D(1), right=NONE),
     dict(alpha=z, beta=lambda x: 1, gamma=lambda x: 2, f=z, x0=0, x1=1, left=('dirichlet', 1), right=('none',), order=1),
     lambda x: exp(-2 * x), lambda x: -2 * exp(-2 * x), [1.8, 1.8, 3.6, 3.6], 'çift p optimal değil'),
    ('T4', '1. mertebe, sağ uç', dict(alpha='', beta='1', gamma='', f='2x', x0=0, x1=2, left=NONE, right=D(4)),
     dict(alpha=z, beta=lambda x: 1, gamma=z, f=lambda x: 2 * x, x0=0, x1=2, left=('none',), right=('dirichlet', 4), order=1),
     lambda x: x ** 2, lambda x: 2 * x, [1.8, 'exact', 'exact', 'exact'], ''),
    ('T5', '2. mertebe BVP', dict(alpha='-1', beta='', gamma='', f='pi^2 sin(pi x)', x0=0, x1=1, left=D(0), right=D(0)),
     dict(alpha=lambda x: -1, beta=z, gamma=z, f=lambda x: pi ** 2 * sin(pi * x), x0=0, x1=1, left=('dirichlet', 0), right=('dirichlet', 0), order=2),
     lambda x: sin(pi * x), lambda x: pi * cos(pi * x), [1.8, 2.8, 3.8, 4.8], ''),
    ('T6', 'karışık D–N', dict(alpha='1', beta='', gamma='-4', f='0', x0=0, x1=1, left=D(0), right=Nm(1)),
     dict(alpha=lambda x: 1, beta=z, gamma=lambda x: -4, f=z, x0=0, x1=1, left=('dirichlet', 0), right=('neumann', 1), order=2),
     lambda x: sinh(2 * x) / (2 * cosh(2)), lambda x: cosh(2 * x) / cosh(2), [1.8, 2.8, 3.8, 4.8], ''),
    ('T7', 'değişken katsayı', dict(alpha='1+x^2', beta='sin(x)', gamma='exp(x)',
                                    f='2x(pi cos(pi x)+2x)+(1+x^2)(-pi^2 sin(pi x)+2)+sin(x)(pi cos(pi x)+2x)+exp(x)(sin(pi x)+x^2)',
                                    x0=0, x1=1, left=D(0), right=D(1)),
     dict(alpha=lambda x: 1 + x ** 2, beta=lambda x: sin(x), gamma=lambda x: exp(x), f=f7, x0=0, x1=1, left=('dirichlet', 0), right=('dirichlet', 1), order=2),
     U7, dU7, [1.8, 2.8, 3.8, 4.8], 'üretilmiş çözüm'),
    ('T8', 'taşınım baskın', dict(alpha='-0.05', beta='1', gamma='', f='1', x0=0, x1=1, left=D(0), right=D(0)),
     dict(alpha=lambda x: -mp.mpf('0.05'), beta=lambda x: 1, gamma=z, f=lambda x: 1, x0=0, x1=1, left=('dirichlet', 0), right=('dirichlet', 0), order=2),
     lambda x: x - (exp(20 * x) - 1) / (exp(20) - 1), lambda x: 1 - 20 * exp(20 * x) / (exp(20) - 1), [None, None, None, None], 'ön-asimptotik; hız iddiası yok'),
    ('T9', 'reaksiyon–difüzyon', dict(alpha='-1', beta='', gamma='100', f='100', x0=0, x1=1, left=D(0), right=D(0)),
     dict(alpha=lambda x: -1, beta=z, gamma=lambda x: 100, f=lambda x: 100, x0=0, x1=1, left=('dirichlet', 0), right=('dirichlet', 0), order=2),
     lambda x: 1 - cosh(10 * (x - mp.mpf(1) / 2)) / cosh(5), lambda x: -10 * sinh(10 * (x - mp.mpf(1) / 2)) / cosh(5), [1.8, 2.8, 3.8, 4.8], ''),
    ('T10', 'Helmholtz tipi, k=1', dict(alpha='1', beta='', gamma='1', f='x', x0=0, x1=1, left=D(0), right=D(0)),
     dict(alpha=lambda x: 1, beta=z, gamma=lambda x: 1, f=lambda x: x, x0=0, x1=1, left=('dirichlet', 0), right=('dirichlet', 0), order=2),
     lambda x: x - sin(x) / sin(1), lambda x: 1 - cos(x) / sin(1), [1.8, 2.8, 3.8, 4.8], ''),
    ('T11', 'iki uç doğal, γ≠0', dict(alpha='1', beta='', gamma='-1', f='-(1+pi^2) cos(pi x)', x0=0, x1=1, left=NONE, right=NONE),
     dict(alpha=lambda x: 1, beta=z, gamma=lambda x: -1, f=lambda x: -(1 + pi ** 2) * cos(pi * x), x0=0, x1=1, left=('none',), right=('none',), order=2),
     lambda x: cos(pi * x), lambda x: -pi * sin(pi * x), [1.8, 2.8, 3.8, 4.8], ''),
    ('T12', 'dejenere uç, serbest', dict(alpha='x', beta='', gamma='', f='1', x0=0, x1=1, left=NONE, right=D(0)),
     dict(alpha=lambda x: x, beta=z, gamma=z, f=lambda x: 1, x0=0, x1=1, left=('none',), right=('dirichlet', 0), order=2),
     lambda x: x - 1, lambda x: 1, ['exact', 'exact', 'exact', 'exact'], 'sınırlı tek çözüm'),
    ('T13', 'dejenere uç, serbest', dict(alpha='x', beta='', gamma='', f='2x', x0=0, x1=1, left=NONE, right=D(0)),
     dict(alpha=lambda x: x, beta=z, gamma=z, f=lambda x: 2 * x, x0=0, x1=1, left=('none',), right=('dirichlet', 0), order=2),
     lambda x: (x ** 2 - 1) / 2, lambda x: x, [1.8, 'exact', 'exact', 'exact'], ''),
    ('T14', '2. mertebe IVP', dict(alpha='1', beta='', gamma='1', f='0', x0=0, x1=1, left=BOTH(0, 1), right=NONE),
     dict(alpha=lambda x: 1, beta=z, gamma=lambda x: 1, f=z, x0=0, x1=1, left=('both', 0, 1), right=('none',), order=2),
     sin, cos, [1.8, 2.8, 3.8, 4.8], 'gözlenen; kanıt yok'),
    ('T15', '2. mertebe IVP, değişken α', dict(alpha='1+x', beta='', gamma='', f='0', x0=0, x1=1, left=BOTH(0, 1), right=NONE),
     dict(alpha=lambda x: 1 + x, beta=z, gamma=z, f=z, x0=0, x1=1, left=('both', 0, 1), right=('none',), order=2),
     lambda x: log(1 + x), lambda x: 1 / (1 + x), [1.8, 2.8, 3.8, 4.8], 'gözlenen; kanıt yok'),
    ('T16', 'Helmholtz k=10', dict(alpha='1', beta='', gamma='100', f='-100', x0=0, x1=1, left=D(0), right=D(0)),
     dict(alpha=lambda x: 1, beta=z, gamma=lambda x: 100, f=lambda x: -100, x0=0, x1=1, left=('dirichlet', 0), right=('dirichlet', 0), order=2),
     lambda x: -1 + cos(10 * x) + Bh * sin(10 * x), lambda x: -10 * sin(10 * x) + 10 * Bh * cos(10 * x), [None, 2.8, 3.8, 4.8], 'P1 ön-asimptotik (kh büyük)'),
]

# Hız iddiasının geçerli olduğu en küçük ağ (asimptotik bölge başlangıcı). Sınır tabakası kalınlığı ~0.1 olan T9'da
# n < 16 ön-asimptotiktir; o ağlarda hız iddia edilmez (tolerans gevşetilmez, iddia kapsamı daraltılır).
ASYMPTOTIC_FROM = {'T9': 16}

# İyi tanımlılık sınır durumları. expected: sürekli sınıf (tek / çözüm yok / tek değil / belirlenemedi) ya da 'ret'.
EDGE = [
    ('E1', dict(alpha='1', beta='', gamma='pi^2', f='sin(pi x)', x0=0, x1=1, left=D(0), right=D(0)), 'no_solution', 'rezonans, ∫f sin πx ≠ 0'),
    ('E1b', dict(alpha='1', beta='', gamma='pi^2', f='0', x0=0, x1=1, left=D(0), right=D(0)), 'non_unique', 'boş uzay sin πx'),
    ('E1c', dict(alpha='1', beta='', gamma='pi^2', f='sin(2 pi x)', x0=0, x1=1, left=D(0), right=D(0)), 'non_unique', 'uyumlu rezonans'),
    ('E1d', dict(alpha='1', beta='', gamma='pi^2+1/1000000', f='sin(pi x)', x0=0, x1=1, left=D(0), right=D(0)), 'unique', 'rezonansa yakın ama tek'),
    ('E2', dict(alpha='x', beta='', gamma='', f='1', x0=0, x1=1, left=D(0), right=D(0)), 'no_solution', 'sınırlı çözüm sınıfında'),
    ('E3', dict(alpha='x-0.5', beta='', gamma='', f='1', x0=0, x1=1, left=D(0), right=D(1)), 'any', 'işaret değiştiren α: uyarı'),
    ('E4', dict(alpha='1', beta='', gamma='', f='2x-1', x0=0, x1=1, left=Nm(0), right=Nm(0)), 'non_unique', 'saf Neumann, uyumlu'),
    ('E4b', dict(alpha='1', beta='', gamma='', f='1', x0=0, x1=1, left=Nm(0), right=Nm(0)), 'no_solution', 'saf Neumann, uyumsuz'),
]
