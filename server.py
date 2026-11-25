#!/usr/bin/env python3
"""server.py — yerel sunucu (yalnızca standart kütüphane).

    python3 server.py            →  http://localhost:8000

Görevi küçük ve net: sayfayı sunar, çözüm kâğıdını fem_sheet.py'ye hesaplatır, istenirse PDF'e derler.
  GET  /            → index.html
  GET  /basis_core.js, /fem_sim.js, /fem_sim3d.js → Basis functions kartı, FEM Simulation (2B) ve 3B görünümü
  POST /api/sheet   → { text, tex, nodal, exact }   (fem_sheet.py ayrı süreçte, süre sınırıyla çalışır)
  POST /api/pdf     → aynı girdiden üretilen .tex'in pdflatex ile derlenmiş hali
  POST /api/exact   → { exact_js }  yalnızca analitik çözüm (grafiklerde bağımsız referans olarak kullanılır)

Neden ayrı süreç: sembolik integral bazı girdilerde çok uzun sürebilir. Süre aşılırsa hesap
"symbolic: false" ile (integraller sayısal) yeniden denenir; sunucu hiçbir girdide kilitlenmez.
Güvenlik: yalnızca 127.0.0.1'i dinler; PDF için istemciden LaTeX metni KABUL EDİLMEZ, .tex her zaman
sunucuda aynı girdiden yeniden üretilir (istemcinin yolladığı keyfi LaTeX derlenmez).
"""
import json
import shutil
import subprocess
import sys
import tempfile
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path

ROOT = Path(__file__).resolve().parent
PORT = 8000
SYMBOLIC_TIMEOUT, NUMERIC_TIMEOUT, LATEX_TIMEOUT, EXACT_TIMEOUT = 60, 60, 40, 20
LATEX = shutil.which('pdflatex')
cache = {}          # girdi (JSON metni) → kâğıt; aynı kâğıt için PDF istenince yeniden hesaplanmasın


def compute_sheet(payload):
    """Kâğıdı hesaplar. Dönüş: (sonuç sözlüğü, HTTP durum kodu)."""
    key = json.dumps(payload, sort_keys=True)
    if key in cache:
        return cache[key], 200
    for symbolic, limit in ((True, SYMBOLIC_TIMEOUT), (False, NUMERIC_TIMEOUT)):
        try:
            run = subprocess.run([sys.executable, str(ROOT / 'fem_sheet.py')], input=json.dumps({**payload, 'symbolic': symbolic}),
                                 capture_output=True, text=True, timeout=limit, encoding='utf-8')
        except subprocess.TimeoutExpired:
            continue                                   # sembolik hesap yetişmedi → sayısal dene
        try:
            result = json.loads(run.stdout)
        except ValueError:
            return {'error': 'The solution sheet could not be generated: ' + (run.stderr.strip().splitlines() or ['unknown error'])[-1]}, 500
        if 'error' in result:
            return result, 400                         # kullanıcının düzeltebileceği girdi hatası
        if len(cache) >= 16:
            cache.pop(next(iter(cache)))
        cache[key] = result
        return result, 200
    return {'error': 'The computation exceeded the time limit. Reduce the number of elements or the polynomial degree.'}, 504


def compute_exact(payload):
    """Analitik çözüm isteği. Bulunamaması hata değildir ({exact_js: None}); grafikler o zaman ince ağ referansına döner."""
    key = json.dumps({**payload, 'mode': 'exact'}, sort_keys=True)
    if key in cache:
        return cache[key]
    try:
        run = subprocess.run([sys.executable, str(ROOT / 'fem_sheet.py')], input=json.dumps({**payload, 'mode': 'exact'}),
                             capture_output=True, text=True, timeout=EXACT_TIMEOUT, encoding='utf-8')
        result = json.loads(run.stdout)
    except (subprocess.TimeoutExpired, ValueError):
        result = {'exact_js': None}
    if 'error' in result:
        result = {'exact_js': None, 'error': result['error']}
    if len(cache) >= 32:
        cache.pop(next(iter(cache)))
    cache[key] = result
    return result


def compile_pdf(tex):
    """LaTeX metnini geçici klasörde derler. Kabuk kullanılmaz, shell-escape kapalıdır."""
    with tempfile.TemporaryDirectory() as folder:
        (Path(folder) / 'sheet.tex').write_text(tex, encoding='utf-8')
        subprocess.run([LATEX, '-no-shell-escape', '-interaction=nonstopmode', '-halt-on-error', 'sheet.tex'],
                       cwd=folder, capture_output=True, timeout=LATEX_TIMEOUT)
        pdf = Path(folder) / 'sheet.pdf'
        return pdf.read_bytes() if pdf.exists() else None


class Handler(BaseHTTPRequestHandler):
    def send(self, status, body, content_type='application/json; charset=utf-8'):
        if isinstance(body, (dict, list)):
            body = json.dumps(body, ensure_ascii=False).encode('utf-8')
        self.send_response(status)
        self.send_header('Content-Type', content_type)
        self.send_header('Content-Length', str(len(body)))
        self.end_headers()
        self.wfile.write(body)

    def do_GET(self):
        if self.path in ('/', '/index.html'):
            self.send(200, (ROOT / 'index.html').read_bytes(), 'text/html; charset=utf-8')
        elif self.path in ('/basis_core.js', '/fem_sim.js', '/fem_sim3d.js'):   # sayfanın eklenti betikleri; sabit beyaz liste, yol girdisi kabul edilmez
            self.send(200, (ROOT / self.path[1:]).read_bytes(), 'text/javascript; charset=utf-8')
        else:
            self.send(404, {'error': 'not found'})

    def do_POST(self):
        try:
            payload = json.loads(self.rfile.read(int(self.headers.get('Content-Length', 0))))
        except ValueError:
            return self.send(400, {'error': 'Invalid request.'})
        if not isinstance(payload, dict):
            return self.send(400, {'error': 'Invalid request.'})
        if self.path == '/api/exact':
            return self.send(200, compute_exact(payload))
        if self.path == '/api/sheet':
            result, status = compute_sheet(payload)
            return self.send(status, result)
        if self.path == '/api/pdf':
            if not LATEX:
                return self.send(501, {'error': 'pdflatex is not installed.'})
            result, status = compute_sheet(payload)
            if status != 200:
                return self.send(status, result)
            try:
                pdf = compile_pdf(result['tex'])
            except subprocess.TimeoutExpired:
                pdf = None
            return self.send(200, pdf, 'application/pdf') if pdf else self.send(500, {'error': 'LaTeX compilation failed.'})
        self.send(404, {'error': 'not found'})

    def log_message(self, *args):       # her istekte terminali doldurmasın
        pass


if __name__ == '__main__':
    print(f'FEM LAB: http://localhost:{PORT}   (PDF: {"pdflatex found" if LATEX else "no pdflatex, print from the browser"})')
    ThreadingHTTPServer(('127.0.0.1', PORT), Handler).serve_forever()
