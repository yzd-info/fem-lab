/* =====================================================================
   fem_sim3d.js — FEM LAB "FEM Simulation" 3B görünümü (Three.js, yalnız istenince yüklenir)

   Neden ayrı ES modülü: Three.js büyüktür ve WebGL gerektirir; çözücü ve 2B görünüm ona hiç bağlı olmasın. index.html bu
   dosyayı yalnız kullanıcı "3D" düğmesine basınca import() ile yükler; yüklenemezse (ağ yok, WebGL yok) 2B'de kalınır.
   Hesap yapmaz: fem_sim.js'in veri sözleşmesini (buildSimulationData) ve 2B ile ORTAK kare durumunu (frameState,
   valueRange) okur. Böylece 2B ve 3B aynı aşamada aynı eleman/terim/değer aralığını gösterir.

   Üçüncü eksen: durağan 1B problemde FİZİKSEL BİR BOYUT DEĞİLDİR. z = −j·Δz, "taban indisi j" katmanıdır: her global
   taban fonksiyonu N_j (ve katkısı a_j N_j) kendi katmanında, toplam u_h ön katmanda (z = 0) durur. Eksen etiketi bunu
   açıkça yazar. Zamana bağlı çözücü eklendiğinde z yerine gerçek t ekseni (temporal.times) kullanılacaktır.

   Hareketli noktalar yalnız anlamlı olanlardır: yaklaşım aşamasında yeni eklenen a_j katsayısı (x_j, a_j) noktası kendi
   katmanından ön katmana (toplama) taşınır. y = a_j sabittir; yalnız z değişir: ara sayısal durum uydurulmaz.
   Çizim isteğe bağlıdır (sürekli döngü yok): kamera değişince, kare değişince ve kısa taşıma animasyonunda çizilir.

   Bölüm haritası: 1. kurulum ve kaynak yönetimi · 2. dünya koordinatları ve malzeme · 3. aşama sahneleri · 4. seçim ·
                   5. dış arayüz (update, resetCamera, dispose)
   ===================================================================== */
import * as THREE from 'three';
import { OrbitControls } from 'three/addons/controls/OrbitControls.js';
import { Line2 } from 'three/addons/lines/Line2.js';
import { LineMaterial } from 'three/addons/lines/LineMaterial.js';
import { LineGeometry } from 'three/addons/lines/LineGeometry.js';

const LX = 10, LY = 4, DEPTH = 6;                       // dünya ölçüleri: x ekseni 10 birim, değer ekseni 4 birim
const canvasHeight = (w) => (w < 500 ? 280 : 360);         // dar ekranda daha alçak: boş dikey alan kalmasın
const MOVE_MS = 450;                                     // a_j noktasının katmandan ön katmana taşınma süresi

/* Neden fabrika: index.html bir kez oluşturur, kare değiştikçe update() çağırır, 2B'ye dönünce dispose() ile GPU
   kaynaklarını bırakır. onPick: düğüm/eleman seçimi durum nesnesine index.html'de yazılır (tek durum nesnesi kuralı). */
export function createSim3D(host, { sim, onPick }) {
  /* ---------- 1. KURULUM VE KAYNAK YÖNETİMİ ---------- */
  const probe = document.createElement('canvas');
  if (!(probe.getContext('webgl2') || probe.getContext('webgl'))) throw new Error('WebGL is not available in this browser');
  const renderer = new THREE.WebGLRenderer({ antialias: true, alpha: true });
  renderer.setPixelRatio(Math.min(window.devicePixelRatio || 1, 2));
  renderer.domElement.setAttribute('role', 'img');
  host.replaceChildren(renderer.domElement);
  const scene = new THREE.Scene(), root = new THREE.Group();
  scene.add(root);
  const camera = new THREE.PerspectiveCamera(32, 1, 0.1, 500);
  const controls = new OrbitControls(camera, renderer.domElement);
  controls.enableDamping = false;                        // damping sürekli çizim döngüsü ister; isteğe bağlı çizimde kapalı
  controls.addEventListener('change', requestRender);
  const lineMaterials = new Set();                       // kalın çizgi malzemeleri: ekran çözünürlüğü değişince güncellenir
  let frameId = 0, moveId = 0, lastKey = '', depth = 1, disposed = false;
  let moving = [];                                       // taşınan a_j noktaları: { mesh, z0, z1 }
  let growing = [];                                      // bu adımda YENİ gelen nesneler: çizgiler uçtan uca büyür, işaretler ölçeklenir
  let lastStepKey = '';

  function requestRender() {
    if (frameId || disposed) return;
    frameId = requestAnimationFrame(() => { frameId = 0; renderer.render(scene, camera); });
  }
  function resize() {
    const w = Math.max(260, host.clientWidth), h = canvasHeight(w);
    renderer.setSize(w, h);
    camera.aspect = w / h;
    camera.updateProjectionMatrix();
    for (const m of lineMaterials) m.resolution.set(w, h);
    requestRender();
  }
  const observer = new ResizeObserver(resize);
  observer.observe(host);

  /* Her kare yeniden kurulur; eski geometri, malzeme ve dokular GPU'dan BIRAKILIR (sızıntı olmasın). */
  function clearScene() {
    cancelAnimationFrame(moveId); moveId = 0; moving = []; growing = [];
    root.traverse((o) => {
      o.geometry?.dispose();
      for (const m of [].concat(o.material || [])) { m.map?.dispose(); m.dispose(); }
    });
    root.clear();
    lineMaterials.clear();
  }

  /* ---------- 2. DÜNYA KOORDİNATLARI VE MALZEME ---------- */
  const css = (name) => getComputedStyle(document.documentElement).getPropertyValue(name).trim() || '#888';
  const color = (name) => new THREE.Color(css(name));
  function fatLine(points, colorName, width, { dashed = false, opacity = 1 } = {}) {
    const flat = [];
    for (const p of points) flat.push(p.x, p.y, p.z);
    const geometry = new LineGeometry();
    geometry.setPositions(flat);
    const material = new LineMaterial({ color: color(colorName), linewidth: width, dashed, dashSize: 0.22, gapSize: 0.14,
      transparent: opacity < 1, opacity, worldUnits: false });
    material.resolution.set(renderer.domElement.width / renderer.getPixelRatio(), renderer.domElement.height / renderer.getPixelRatio());
    lineMaterials.add(material);
    const line = new Line2(geometry, material);
    if (dashed) line.computeLineDistances();
    return line;
  }
  function thinLine(points, colorName, opacity = 1) {
    const geometry = new THREE.BufferGeometry().setFromPoints(points);
    return new THREE.Line(geometry, new THREE.LineBasicMaterial({ color: color(colorName), transparent: opacity < 1, opacity }));
  }
  function sphere(pos, r, colorName, { hollow = false, userData = null } = {}) {
    const mesh = new THREE.Mesh(new THREE.SphereGeometry(r, 16, 12),
      new THREE.MeshBasicMaterial({ color: color(hollow ? '--surface' : colorName) }));
    mesh.position.copy(pos);
    if (hollow) {                                        // boş düğüm (eleman içi): her açıdan görünen renkli dış çeper (arka yüz)
      mesh.add(new THREE.Mesh(new THREE.SphereGeometry(r * 1.4, 16, 12), new THREE.MeshBasicMaterial({ color: color(colorName), side: THREE.BackSide })));
    }
    if (userData) mesh.userData = { ...mesh.userData, ...userData };
    return mesh;
  }
  /* Yazı: tuval dokulu sprite (ekrana dönük). Sayı az tutulur (seyreltme), her biri dispose edilir. */
  function label(text, pos, { colorName = '--ink-2', size = 0.34, align = 'center' } = {}) {
    const c = document.createElement('canvas'), ctx = c.getContext('2d'), px = 48;
    ctx.font = `${px}px "IBM Plex Mono", ui-monospace, monospace`;
    c.width = Math.ceil(ctx.measureText(text).width) + 16; c.height = px + 16;
    ctx.font = `${px}px "IBM Plex Mono", ui-monospace, monospace`;
    ctx.fillStyle = css(colorName); ctx.textBaseline = 'middle';
    ctx.fillText(text, 8, c.height / 2);
    const texture = new THREE.CanvasTexture(c);
    texture.colorSpace = THREE.SRGBColorSpace;
    const sprite = new THREE.Sprite(new THREE.SpriteMaterial({ map: texture, transparent: true, depthWrite: false }));
    const h = size, w = size * c.width / c.height;
    sprite.scale.set(w, h, 1);
    sprite.center.set(align === 'left' ? 0 : align === 'right' ? 1 : 0.5, 0.5);
    sprite.position.copy(pos);
    return sprite;
  }
  function box(xa, xb, ya, yb, za, zb, colorName, opacity, userData = null) {
    const mesh = new THREE.Mesh(new THREE.BoxGeometry(xb - xa, yb - ya, Math.max(zb - za, 0.02)),
      new THREE.MeshBasicMaterial({ color: color(colorName), transparent: true, opacity, depthWrite: false }));
    mesh.position.set((xa + xb) / 2, (ya + yb) / 2, (za + zb) / 2);
    if (userData) { mesh.userData = userData; mesh.material.visible = opacity > 0; }
    return mesh;
  }

  /* ---------- 3. AŞAMA SAHNELERİ ---------- */
  function build(data, stages, step, opts) {
    const fs = sim.frameState(data, stages, step, opts), { stage, local: k, selNode, selEl } = fs;
    const D = data.discretization, S = data.solution, B = data.basis;
    const { elementCount: n, nodeCount, polynomialDegree: p, coordinates: X, connectivity: C, elementBoundaries: EB } = D;
    const [x0, x1] = data.problem.domain;
    const range = sim.valueRange(data, stage.id);
    const layered = ['local', 'global', 'approx'].includes(stage.id);
    depth = layered ? DEPTH : 1;
    const dz = layered ? DEPTH / Math.max(1, nodeCount) : 0;
    const wx = (x) => (x - x0) / (x1 - x0) * LX - LX / 2;
    const wy = range ? (v) => (v - range.lo) / (range.hi - range.lo) * LY : () => LY * 0.35;
    const zLayer = (j) => -j * dz;
    const y0 = range ? wy(0) : LY * 0.35;
    const V = (x, y, z = 0) => new THREE.Vector3(wx(x), y, z);
    const curve = (xs, vs, z) => {                        // eksik değer çizgiyi böler (0 sanılmaz)
      const parts = [[]];
      xs.forEach((x, i) => { const v = vs[i]; if (v === null || !Number.isFinite(v)) { if (parts.at(-1).length) parts.push([]); } else parts.at(-1).push(V(x, wy(v), z)); });
      return parts.filter((q) => q.length > 1);
    };
    const add = (o) => { root.add(o); return o; };
    const grow = (o) => { growing.push(o); return o; };   // yeni gelen nesne (2B'deki sweep/enter karşılığı)
    const counts = { nodes: 0, layers: 0, elements: 0, moving: 0 };

    /* Eksenler: x (uzay), değer (aşamanın kendi ölçeği), katmanlı aşamada z = taban indisi j (fiziksel değil). */
    add(thinLine([V(x0, y0), V(x1, y0)], '--c-mesh'));
    if (opts.showLabels !== false) {
      for (const t of sim.xTicks(x0, x1, 4, !!data.problem.angularX)) add(label(t.label, V(t.value, y0 - 0.32, 0.25), { size: 0.3 }));
      add(label(data.problem.angularX ? 'x (rad) →' : 'x →', V(x1, y0 - 0.32, 0.25).add(new THREE.Vector3(0.45, 0, 0)), { size: 0.3, align: 'left' }));
    }
    if (range) {
      add(thinLine([V(x0, 0), V(x0, LY)], '--c-mesh'));
      if (opts.showLabels !== false) {
        for (const t of sim.niceTicks(range.lo, range.hi, 4)) add(label(sim.fmt(t), V(x0, wy(t)).add(new THREE.Vector3(-0.2, 0, 0)), { size: 0.28, align: 'right' }));
        const title = { local: 'N_k(x)', global: 'N_j(x)', approx: 'a_j N_j(x), u_h(x)', reference: 'u(x)', error: data.metadata.errorKind === 'exact' ? 'e = u − u_h' : 'd = u_ref − u_h' }[stage.id];
        add(label(`↑ ${title}`, V(x0, LY + 0.35), { size: 0.28, align: 'left', colorName: '--ink' }));
      }
    }
    if (layered) {
      add(thinLine([V(x0, y0, 0), V(x0, y0, zLayer(nodeCount))], '--c-mesh'));
      if (opts.showLabels !== false) add(label('depth: basis index j (not a physical dimension)', V(x0, y0 + 0.25, zLayer(nodeCount) - 0.3), { size: 0.3, align: 'left' }));
    }
    /* Ağ: eleman sınırları ön katmanda ince dikey çizgiler; seçim için görünmez eleman kutuları. */
    if (opts.showMesh !== false && stage.id !== 'domain') for (let e = 0; e <= fs.shownElements; e++)
      add(thinLine([V(EB[e], range ? 0 : y0 - 0.25), V(EB[e], range ? LY : y0 + 0.25)], '--c-mesh', 0.45));
    if (stage.id !== 'domain') for (let e = 1; e <= fs.shownElements; e++)
      add(box(wx(EB[e - 1]), wx(EB[e]), range ? 0 : y0 - 0.4, range ? LY : y0 + 0.4, -0.05, 0.05, '--c-mesh', 0, { element: e }));
    const nodeSphere = (j, y, z, colorName, r = 0.07) => {
      counts.nodes++;
      return add(sphere(V(X[j - 1], y, z), r, colorName, { hollow: (j - 1) % p !== 0, userData: { node: j } }));
    };
    const labelEvery = (count) => Math.max(1, Math.ceil(count / 12));   // en çok ~12 numara: sprite sayısı ve okunurluk

    if (stage.id === 'domain') {
      const line = add(fatLine([V(x0, y0), V(x1, y0)], '--ink', 4));
      if (k === 1) grow(line);
      if (opts.showLabels !== false) {
        add(label(`x₀ = ${sim.fmt(x0)}`, V(x0, y0 - 0.45), { align: 'left' }));
        add(label(`x₁ = ${sim.fmt(x1)}`, V(x1, y0 - 0.45), { align: 'right' }));
        if (k >= 2) {
          const [L, R] = data.problem.boundaryConditions;
          const l1 = add(label(sim.conditionLabel(L, data.problem), V(x0, y0 + 0.6), { align: 'left', colorName: '--ink' }));
          const l2 = add(label(sim.conditionLabel(R, data.problem), V(x1, y0 + 0.6), { align: 'right', colorName: '--ink' }));
          if (k === 2) { grow(l1); grow(l2); }
        }
      }
      if (k === 0) {
        /* Girdi fonksiyonları (ilk kare): her biri kendi katmanında ve KENDİ ölçeğinde (2B'deki küçük çoklu grafiklerin karşılığı). */
        (data.problem.inputFunctions || []).forEach((fn, i) => {
          const z = -(i + 1) * 1.8;
          let lo = Infinity, hi = -Infinity;
          for (const v of fn.values) if (v !== null) { lo = Math.min(lo, v); hi = Math.max(hi, v); }
          if (!Number.isFinite(lo)) { lo = -1; hi = 1; }
          if (hi - lo < 1e-12 * Math.max(1, Math.abs(hi))) { const c = Math.max(Math.abs(hi), 1) * 0.5; lo -= c; hi += c; }
          const fy = (v) => (v - lo) / (hi - lo) * LY * 0.8;
          add(thinLine([V(x0, lo <= 0 && hi >= 0 ? fy(0) : 0, z), V(x1, lo <= 0 && hi >= 0 ? fy(0) : 0, z)], '--c-mesh', 0.6));
          const pts = [];
          S.samples.x.forEach((x, s2) => { const v = fn.values[s2]; if (v !== null) pts.push(V(x, fy(v), z)); });
          if (pts.length > 1) grow(add(fatLine(pts, '--c-input', 3)));
          /* Etiket eğrinin sağ ucunda, kendi katmanında: katmanlar üst üste binse de etiketler ayrışır. */
          if (opts.showLabels !== false && pts.length) add(label(fn.constant !== null ? `${fn.name} = ${sim.fmt(fn.constant)}` : fn.name, pts.at(-1).clone().add(new THREE.Vector3(0.2, 0.2, 0)), { align: 'left', colorName: '--c-input' }));
        });
      }
    } else if (stage.id === 'mesh' || stage.id === 'element') {
      add(fatLine([V(x0, y0), V(EB[fs.shownElements], y0)], '--ink', 3));
      counts.elements = fs.shownElements;
      if (stage.id === 'element') { const b = add(box(wx(EB[selEl - 1]), wx(EB[selEl]), y0 - 0.5, y0 + 0.5, -0.3, 0.3, '--c-fem', 0.12)); if (k === 0) grow(b); }
      if (opts.showLabels !== false) {                    // eleman numaraları (2B ile aynı bilgi), en çok ~12
        const everyEl = labelEvery(fs.shownElements);
        for (let e = 1; e <= fs.shownElements; e++) if ((e - 1) % everyEl === 0 || e === fs.shownElements || (stage.id === 'element' && e === selEl))
          add(label(`e${e}`, V((EB[e - 1] + EB[e]) / 2, y0 - 0.38), { size: 0.26, colorName: stage.id === 'element' && e === selEl ? '--c-fem' : '--c-mesh' }));
      }
      if (stage.id === 'element' && k >= 1) {
        /* Referans eleman s ∈ [0, h] (kâğıttaki gösterim) eksenin altında; kesikli çizgiler yerel düğüm k → global düğüm (e − 1)p + k. */
        const yRef = y0 - 1.5, mid = (wx(EB[selEl - 1]) + wx(EB[selEl])) / 2, half = Math.max((wx(EB[selEl]) - wx(EB[selEl - 1])) / 2, 1.4);
        const ra = Math.max(-LX / 2, Math.min(mid - half, LX / 2 - 2 * half)), rx = (t) => ra + (t + 1) * half;
        grow(add(fatLine([new THREE.Vector3(ra, yRef, 0), new THREE.Vector3(ra + 2 * half, yRef, 0)], '--ink-2', 2)));
        if (opts.showLabels !== false) add(label('s', new THREE.Vector3(ra - 0.25, yRef, 0), { align: 'right' }));
        B.referenceNodes.forEach((t, kk) => {
          const j = C[selEl - 1][kk];
          grow(add(thinLine([new THREE.Vector3(rx(t), yRef, 0), V(X[j - 1], y0, 0)], '--c-mesh', 0.7)));
          grow(add(sphere(new THREE.Vector3(rx(t), yRef, 0), 0.06, '--ink-2', { hollow: kk !== 0 && kk !== p })));
          if (opts.showLabels !== false) {
            add(label(sim.fmt(D.h * (1 + t) / 2), new THREE.Vector3(rx(t), yRef - 0.32, 0), { size: 0.24 }));
            add(label(`${kk + 1}→${j}`, new THREE.Vector3(rx(t), yRef + 0.3, 0), { size: 0.24 }));
          }
        });
      }
      const every = labelEvery(nodeCount);
      for (let j = 1; j <= nodeCount; j++) {
        const end = (j - 1) % p === 0;
        if (X[j - 1] > EB[fs.shownElements] + 1e-12 * (x1 - x0) || (!end && !fs.finalMesh)) continue;
        const inSel = stage.id === 'element' && C[selEl - 1].includes(j);
        const isNew = stage.id === 'mesh' && X[j - 1] > EB[fs.prevShown] + 1e-12 * (x1 - x0), lastMesh = stage.id === 'mesh' && k === stage.steps - 1;
        const dot = nodeSphere(j, y0, 0, inSel ? '--c-fem' : '--ink', inSel ? 0.1 : 0.07);
        if (isNew || (lastMesh && !end)) grow(dot);
        if (fs.finalMesh && opts.showLabels !== false && ((j - 1) % every === 0 || j === nodeCount || inSel)) {
          const lab = add(label(String(j), V(X[j - 1], y0 + 0.32), { size: 0.3, colorName: inSel ? '--c-fem' : '--ink-2' }));
          if (lastMesh) grow(lab);
        }
      }
    } else if (stage.id === 'local' || stage.id === 'global') {
      /* Her yerel N_k, ait olduğu global düğümün katmanında: katmandaki parçaların birleşimi N_j'dir. */
      const piece = (e, kk, colorName, width, opacity = 1, fresh = false) => {
        const j = C[e - 1][kk];
        for (const pts of curve(B.elementX[e - 1], B.N[kk], zLayer(j))) {
          const o = add(width > 1.2 ? fatLine(pts, colorName, width, { opacity }) : thinLine(pts, colorName, opacity));
          if (fresh) grow(o);
        }
      };
      if (stage.id === 'local') {
        add(box(wx(EB[selEl - 1]), wx(EB[selEl]), 0, LY, zLayer(C[selEl - 1][p]) - dz / 2, 0.05, '--c-basis', 0.06));
        for (let kk = 0; kk <= k; kk++) piece(selEl, kk, kk === k ? '--c-basis' : '--c-basis-soft', kk === k ? 3 : 2, 1, kk === k);
        C[selEl - 1].forEach((j, kk) => nodeSphere(j, wy(kk === k ? 1 : 0), zLayer(j), '--c-basis'));
        counts.layers = k + 1;
      } else {
        const sup = B.support[selNode - 1];
        const near = new Set([selNode, ...fs.neighbors]);
        for (let e = 1; e <= n; e++) for (let kk = 0; kk <= p; kk++) {
          const j = C[e - 1][kk];
          if (k >= 1 && j === selNode) continue;
          piece(e, kk, k >= 2 && near.has(j) ? '--c-basis' : '--c-basis-soft', k >= 2 && near.has(j) ? 1.8 : 1, k >= 2 && near.has(j) ? 0.8 : 1, k === 0 || (k === 2 && near.has(j)));
        }
        if (k >= 1) {
          for (const { element, local } of sup) piece(element, local - 1, '--c-basis', 3.2, 1, k === 1);
          const supBox = add(box(wx(EB[sup[0].element - 1]), wx(EB[sup.at(-1).element]), 0, LY, zLayer(selNode) - dz / 2, zLayer(selNode) + dz / 2, '--c-basis', 0.1));
          if (k === 1) grow(supBox);
        }
        for (let j = 1; j <= nodeCount; j++) if (nodeCount <= 129 || near.has(j)) nodeSphere(j, wy(j === selNode && k >= 1 ? 1 : 0), zLayer(j), j === selNode ? '--c-basis' : '--c-mesh', j === selNode ? 0.1 : 0.05);
        counts.layers = nodeCount;
      }
    } else if (stage.id === 'approx') {
      const a = S.nodalCoefficients;
      /* Katkılar a_j N_j kendi katmanında; kısmi toplam (son adımda u_h) ön katmanda. Kısmi toplam, 2B ile aynı formül. */
      const from = Math.max(1, fs.terms - (nodeCount > 65 ? 64 : nodeCount) + 1);
      for (let j = from; j <= fs.terms; j++) for (const { element, local } of B.support[j - 1]) {
        const vals = B.N[local - 1].map((v) => a[j - 1] * v), fresh = fs.newestFrom !== null && j >= fs.newestFrom;
        for (const pts of curve(B.elementX[element - 1], vals, zLayer(j))) {
          const o = add(fresh || j === selNode ? fatLine(pts, '--c-basis', 2, { dashed: true }) : thinLine(pts, '--c-basis-soft'));
          if (fresh) grow(o);
        }
      }
      counts.layers = fs.terms;
      for (let e = 1; e <= n; e++) {
        const vals = B.elementX[e - 1].map((_, s) => C[e - 1].reduce((acc, j, kk) => acc + (j <= fs.terms ? a[j - 1] * B.N[kk][s] : 0), 0));
        for (const pts of curve(B.elementX[e - 1], vals, 0)) { const o = add(fatLine(pts, '--c-fem', fs.complete ? 4 : 3)); if (!fs.complete) grow(o); }
      }
      const animate = fs.newestFrom !== null && !matchMedia('(prefers-reduced-motion: reduce)').matches;
      for (let j = 1; j <= fs.terms; j++) {
        if (nodeCount > 129 && j !== selNode && !(fs.newestFrom !== null && j >= fs.newestFrom)) continue;
        const fresh = animate && j >= fs.newestFrom;
        const s = nodeSphere(j, wy(a[j - 1]), fresh ? zLayer(j) : 0, '--c-fem', j === selNode ? 0.1 : 0.07);
        if (fresh) { moving.push({ mesh: s, z0: zLayer(j), z1: 0 }); counts.moving++; }
      }
    } else if (stage.id === 'reference') {
      const fem = curve(S.samples.x, S.samples.fem, 0), ref = curve(S.samples.x, S.samples.reference, 0.02);
      for (const pts of fem) add(fatLine(pts, '--c-fem', 4, { opacity: k === 0 ? 0.35 : 1 }));
      for (const pts of ref) { const o = add(fatLine(pts, '--c-ref', 3, { dashed: true })); if (k === 0) grow(o); }
      if (k >= 1 && nodeCount <= 129) for (let j = 1; j <= nodeCount; j++) nodeSphere(j, wy(S.nodalCoefficients[j - 1]), 0, '--c-fem', 0.06);
    } else if (stage.id === 'error') {
      const err = S.samples.error;
      for (const pts of curve(S.samples.x, err, 0)) { const o = add(fatLine(pts, '--c-err', 3)); if (k === 0) grow(o); }
      let peak = 0;
      for (const v of err) if (v !== null) peak = Math.max(peak, Math.abs(v));
      if (k >= 1 && peak > 0) for (let e = 1; e <= n; e++) {   // 2B ile aynı ölçüt: eleman içi en büyük |e| ≥ genel en büyüğün yarısı
        let m = 0;
        S.samples.x.forEach((x, i) => { if (x >= EB[e - 1] - 1e-15 && x <= EB[e] + 1e-15 && err[i] !== null) m = Math.max(m, Math.abs(err[i])); });
        if (m >= 0.5 * peak) { grow(add(box(wx(EB[e - 1]), wx(EB[e]), 0, LY, -0.2, 0.2, '--c-err', 0.12))); counts.elements++; }
      }
    }
    host.dataset.stage = stage.id;                         // test ve ekran okuyucu için: hangi aşama, kaç nesne
    host.dataset.nodes = String(counts.nodes);
    host.dataset.layers = String(counts.layers);
    host.dataset.elements = String(counts.elements);
    host.dataset.moving = String(counts.moving);
    host.dataset.growing = String(growing.length);
    renderer.domElement.setAttribute('aria-label', `3D view, ${stage.title}`);
    return layered;
  }

  /* Büyüme: çizgi yalnız ilk t kesri kadar çizilir (drawRange / instanceCount), işaret ve kutu ölçeklenir. Değerler değişmez. */
  function applyGrow(t) {
    for (const o of growing) {
      if (o.isLine2) { const segs = o.geometry.attributes.instanceStart.count; o.geometry.instanceCount = Math.ceil(segs * t); }
      else if (o.isLine) o.geometry.setDrawRange(0, Math.ceil(o.geometry.attributes.position.count * t));
      else { o.userData.s0 ||= o.scale.clone(); o.scale.copy(o.userData.s0).multiplyScalar(Math.max(t, 1e-3)); }
    }
  }
  /* a_j noktalarını katmandan ön katmana taşır (yalnız z) ve yeni nesneleri büyütür; bittiğinde döngü durur. */
  function runMove(durationMs) {
    if (!moving.length && !growing.length) return;
    const t0 = performance.now(), dur = durationMs || MOVE_MS;
    const tick = (now) => {
      const t = Math.min(1, (now - t0) / dur), ease = t * t * (3 - 2 * t);
      for (const m of moving) m.mesh.position.z = m.z0 + (m.z1 - m.z0) * ease;
      applyGrow(ease);
      renderer.render(scene, camera);
      moveId = t < 1 ? requestAnimationFrame(tick) : 0;
      if (t >= 1) { moving = []; growing = []; }
    };
    moveId = requestAnimationFrame(tick);
  }

  /* ---------- 4. SEÇİM ----------
     Sürükleme (kamera döndürme) seçim sayılmaz: basma ile bırakma arası 5 pikselden azsa tıklamadır. */
  const raycaster = new THREE.Raycaster(), pointer = new THREE.Vector2();
  let down = null;
  renderer.domElement.addEventListener('pointerdown', (e) => { down = [e.clientX, e.clientY]; });
  renderer.domElement.addEventListener('pointerup', (e) => {
    if (!down || Math.hypot(e.clientX - down[0], e.clientY - down[1]) > 5) return;
    const r = renderer.domElement.getBoundingClientRect();
    pointer.set(((e.clientX - r.left) / r.width) * 2 - 1, -((e.clientY - r.top) / r.height) * 2 + 1);
    raycaster.setFromCamera(pointer, camera);
    const hits = raycaster.intersectObjects(root.children, true).map((h) => h.object.userData.node ? h.object : h.object.parent?.userData.node ? h.object.parent : h.object);
    const node = hits.find((o) => o.userData.node), element = hits.find((o) => o.userData.element);
    if (node) onPick({ node: node.userData.node }); else if (element) onPick({ element: element.userData.element });
  });

  /* ---------- 5. DIŞ ARAYÜZ ---------- */
  function resetCamera() {
    /* Önden-soldan-üstten: x ekseni soldan sağa okunur (2B ile aynı yön), katmanlar arkaya doğru görünür. */
    /* Dar ekranda (en-boy < 1.6) yatay görüş açısı daralır: kamera aynı doğrultuda geri çekilir ki x ekseni sığsın. */
    controls.target.set(LX * 0.04, LY * 0.4, -depth / 2);
    const back = Math.max(1, 1.6 / camera.aspect);
    camera.position.set(-LX * 0.28, LY * 1.55, LX * 1.05 + depth * 0.4).sub(controls.target).multiplyScalar(back).add(controls.target);
    controls.update();
    requestRender();
  }
  let cameraSet = false;
  return {
    update(data, stages, step, opts) {
      const key = JSON.stringify([data.version, data.discretization.nodeCount, data.solution.nodalCoefficients.length, step, opts, opts.dataKey, document.documentElement.dataset.theme || '', matchMedia('(prefers-color-scheme: dark)').matches]);
      if (key === lastKey) return;
      lastKey = key;
      clearScene();
      build(data, stages, step, opts);
      if (!cameraSet) { resize(); resetCamera(); cameraSet = true; }
      /* Yalnız adım ya da veri değişince canlandır (etiket/ağ düğmesi aynı kareyi yeniden kurar, animasyon tekrar etmesin). */
      const stepKey = JSON.stringify([opts.dataKey, step]), animate = stepKey !== lastStepKey && !matchMedia('(prefers-reduced-motion: reduce)').matches;
      lastStepKey = stepKey;
      if (!animate) { for (const m of moving) m.mesh.position.z = m.z1; moving = []; growing = []; }
      else applyGrow(0);
      requestRender();
      runMove(opts.durationMs);
    },
    resetCamera() { controls.target.z = -depth / 2; resetCamera(); },
    dispose() {
      disposed = true;
      cancelAnimationFrame(frameId); cancelAnimationFrame(moveId);
      observer.disconnect();
      clearScene();
      controls.dispose();
      renderer.dispose();
      renderer.domElement.remove();
      for (const k of ['stage', 'nodes', 'layers', 'elements', 'moving', 'growing']) delete host.dataset[k];
    }
  };
}
