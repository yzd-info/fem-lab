"""ref_parser.py — doğrulama araçlarının ortak girdi okuyucusu (yalnız SymPy; FEM Lab'ın ayrıştırıcısını kullanmaz).

Neden ayrı dosya: skfem_reference (scikit-fem referansı) ve continuous_problem (sürekli düzey doğrulaması) AYNI okuyucuyu
kullanmalı (üçüncü bir ayrıştırıcı yazılmasın); ama continuous_problem scikit-fem'e bağımlı olmamalı — tests.py onu
çalıştırır ve yalnız requirements.txt kurulu bir makinede de geçmelidir.
"""
import sympy as sp
from sympy.parsing.sympy_parser import parse_expr, standard_transformations, implicit_multiplication_application, convert_xor

X = sp.Symbol('x', real=True)
NAMES = {'x': X, 'pi': sp.pi, 'e': sp.E, 'sin': sp.sin, 'cos': sp.cos, 'tan': sp.tan, 'asin': sp.asin, 'acos': sp.acos,
         'atan': sp.atan, 'sinh': sp.sinh, 'cosh': sp.cosh, 'tanh': sp.tanh, 'exp': sp.exp, 'ln': sp.log, 'log': sp.log,
         'sqrt': sp.sqrt, 'abs': sp.Abs}
GRAMMAR = standard_transformations + (implicit_multiplication_application, convert_xor)


# Girdi metnini SymPy'nin kendi dilbilgisiyle okur (FEM Lab'ın ayrıştırıcısı kullanılmaz). Yalnız FEM Lab'ın KABUL ettiği
# girdilerde çağrılır; güvenlik denetimi FEM Lab'ın işidir (bu modül sunucuda çalışmaz).
def expr(text):
    text = (str(text) if text is not None else '').strip().lower() or '0'
    e = sp.sympify(parse_expr(text, local_dict=dict(NAMES), transformations=GRAMMAR))
    # Ondalık girdi tam kesirdir (uygulamayla aynı anlam: 9.8 = 49/5); Float'lar dsolve'u kararsızlaştırabiliyor.
    return e.xreplace({v: sp.Rational(str(v)) for v in e.atoms(sp.Float)})
