# Architecture

This document describes how FEM LAB is built: the files, the two computation kernels, the five-layer pipeline of the
server, the user interface, and how to extend the application. The mathematics is in [MATEMATIK.md](MATEMATIK.md);
verification and testing in [DOGRULAMA.md](DOGRULAMA.md).

## 1. Files

| File | Role |
|---|---|
| `index.html` | User interface and the **numerical** FEM kernel for the live charts (JavaScript, no dependencies; MathJax is loaded from a CDN) |
| `basis_core.js` | Data for the basis-function views, built on the solver's own shape functions |
| `fem_sim.js`, `fem_sim3d.js` | Simulation views in 2D (DOM-free) and 3D (Three.js, loaded from a CDN); they consume solver data and do not solve |
| `fem_sheet.py` | Layers 1, 2, 3 and 5 of the server pipeline: problem definition, problem validation, the **symbolic/exact** FEM engine (SymPy), and the solution sheet (formatting only) |
| `fem_verify.py` | Layer 4: independent verification; does not call the FEM engine |
| `server.py` | Local web server (Python standard library only): serves the page, runs the sheet computation, compiles PDFs |
| `tests.js` | Tests of the browser kernel (`node tests.js`) |
| `tests_basis.js`, `tests_sim.js` | Tests of the basis-function data and the simulation views |
| `tests.py` | Tests of the server pipeline and the sheet (`python3 tests.py`) |
| `requirements.txt` | Python dependency: `sympy` (which brings `mpmath`) |
| `verification/` | `ifem.py` (an independent mpmath FEM that imports nothing from the application), `benchmarks.py` (T1–T16, E1–E4b), `run_verification.py` (`--quick`, `--full`, `--baseline`, `--compare`), `baseline/` (recorded reference values); `corpus/cases.json` (80 input-only cases), `run_corpus.py` (corpus runner, cross-verification matrix), `skfem_reference.py` (scikit-fem reference), `continuous_problem.py` (continuous existence/uniqueness/compatibility), `ref_parser.py` (their shared SymPy input reader), `corpus_results/` (generated reports; `baseline_v1/` holds the data of an earlier run used for comparison) |
| `.github/workflows/ci.yml` | Continuous integration |
| `ornekler/` | Example solution sheets (PDF) |
| `docs/` | Technical documentation and screenshots |

## 2. Two kernels

The JavaScript kernel solves numerically so that the charts update instantly while sliders are dragged. The SymPy
kernel solves exactly so that every integral on the solution sheet can be written out. The two are deliberately
independent and are compared on every sheet (see [DOGRULAMA.md §1](DOGRULAMA.md#1-two-independent-kernels)).

The page works without the server: the charts use the browser kernel and the fine-mesh reference. The solution sheet,
the PDF export and the analytical reference require the server.

### Server endpoints

| Endpoint | Returns |
|---|---|
| `GET /` | `index.html` |
| `GET /basis_core.js`, `/fem_sim.js`, `/fem_sim3d.js` | the page's add-on scripts (fixed whitelist) |
| `POST /api/sheet` | the solution sheet (text, LaTeX source, nodal values, analytical solution, verification items, …) |
| `POST /api/pdf` | the sheet compiled with `pdflatex` |
| `POST /api/exact` | only the analytical solution and continuous classification (used as an independent reference by the charts) |

The sheet is computed in a separate process with a time limit. If the symbolic computation exceeds 60 s, it is
repeated with numerical integrals, so the server never hangs on an input. The server listens only on `127.0.0.1`, and
for the PDF it never accepts LaTeX from the client: the `.tex` file is always regenerated on the server from the same
input. Results are cached in memory per input; restart the server after changing `fem_sheet.py` or `fem_verify.py`.

## 3. Five-layer pipeline (server)

Data flows only downwards; each layer receives the structured output of the previous one. If a layer fails, the later
layers do not run and the reason is reported (no sheet is produced and no unsupported claim is written).

| Layer | Where | Responsibility | Does not |
|---|---|---|---|
| 1 Problem definition | `fem_sheet.define_problem` | stores the user's text exactly as entered (single source of truth) | interpret or correct the input |
| 2 Problem validation | `fem_sheet.validate_problem` | parses every field with two independent parsers (token translator and SymPy grammar), records the equivalence in a trace table, records the meaning of each boundary condition (essential / derivative / flux / initial) and its flux rule, stops on rule violations; the browser compares a third parser at the same 9 points | ask whether a solution exists or is unique; silently change the input |
| 3 FEM engine | `fem_sheet.fem_engine` | shape-function construction and checks, Galerkin formulation, element integrals, assembly, boundary conditions, solve, reconstruction | produce LaTeX text |
| 4 Independent verification | `fem_verify.verify` | checks the engine's output by independent routes (see [DOGRULAMA.md §6](DOGRULAMA.md#6-independent-verification-layer-fem_verifypy)) | trust the engine's success status |
| 5 Solution sheet | `fem_sheet.build_sheet` | formats the data of layers 2–4 | compute any number |

`run_pipeline` applies this order. If the engine reports a singular discrete system, the continuous classification of
layer 4 is still run: if the continuous problem is symbolically proven to have no or no unique solution, that more
fundamental reason is reported; otherwise the discrete diagnosis is reported together with the continuous information.
The layer boundaries (layer 4 does not call the engine; layer 5 performs no computation) are checked in `tests.py`
from the source code.

## 4. User interface

- **Input panel:** the four coefficients $\alpha, \beta, \gamma, f$, the interval $[x_0, x_1]$, the boundary condition at
  each end, the method, and sliders for $n$ (1–64) and $p$ (1–4). The page starts empty; there are no preset problems.
  While the $n$ slider has focus, the arrow keys and Home/End work; the labels 1 / 8 / 16 / 32 / 64 below it jump to
  that value.
- **Appearance:** the switch at the top right selects Auto (follow the system), Light or Dark; the choice is stored in
  the browser and applied before the first paint. Printing is always black on white.
- **FEM simulation · Approximation:** the continuous problem, mesh generation, element structure, local shape functions
  and global basis construction, step by step.
- **FEM simulation · Solution:** the sections of the solution sheet in order, from the approximate solution through the
  element matrices (cell by cell, with the integrand), assembly, boundary conditions and solve to the error measures.
- Both simulation panels show formulas taken from the solution sheet. When a panel is visible, the page requests the
  sheet from the server (about 0.9 s after the last change). A sheet that does not match the current input is never
  shown; without the server, or for a step the sheet does not contain, a built-in formula is shown. The steps are
  construction steps, not physical time.
- **Summary cards:** mesh, reference type, observed convergence orders, and the $L^2$, $H^1$, maximum and nodal
  measures; clicking a card shows its definition.
- **Basis functions:** global basis functions with a node slider; the details show the derivative, the approximation
  built from the computed coefficients and a derivation.
- **Solution chart:** the reference (blue, long dashes) and $u_h$ for the selected degrees, each evaluated with its own
  shape functions on its own elements (no smoothing). Chips select the degrees. Each degree shows its own nodes
  ($N = np + 1$): filled markers are element end nodes, hollow markers interior nodes. Dashed vertical lines are element
  boundaries. If a coefficient is trigonometric (without `pi` in its argument), the $x$ axis is labelled in multiples
  of $\pi$ when the interval allows it.
- **Pointwise error:** "signed" view $u - u_h$ (linear axis) or "absolute (log)" view $|u - u_h|$ with the nodal errors
  as points. The shaded region is the measurement floor.
- **Convergence (h / p):** "h" is h-refinement (degree fixed, $n$ varies; log–log axes, measured rate at the end of each
  line, theoretical slope as a dashed grey line only for second-order boundary value problems); "p" is p-refinement
  ($n$ fixed, degree 1 → 4; measured values only, no rate because $h$ does not change). hp-refinement is not shown.
  Filled markers contribute to the rate, hollow markers lie at the measurement floor; an undefined rate is shown as
  "—", never as 0. $L^2$ / $H^1$ can be switched. The error and convergence charts can be zoomed and panned.
- **Numerical diagnostics:** continuous classification, degenerate coefficients, first-order compatibility,
  conditioning, quadrature-error estimate and stability indicators. Warnings stay visible; the metrics expand on demand.
- **Solution sheet:** opened and closed by clicking its heading or the "Open solution sheet" button. While it is open,
  a change to the problem, $n$ or $p$ recomputes it after about 0.7 s (at most one request from this panel at a time).
  "Download PDF" returns the compiled PDF (if `pdflatex` is not installed, the browser's print dialog opens instead);
  "Download .tex" and "Copy LaTeX source" give the source. The cross-check and the verification items are listed above
  the sheet.

## 5. Adding a method

1. Add an entry to the `METHODS` object in `index.html` (inside the `/*CORE-START*/ … /*CORE-END*/` block, which must not
   touch the DOM because `tests.js` extracts and runs it):

   ```js
   petrov: {
     label: 'Petrov–Galerkin (SUPG)',
     solve(problem, mesh) {
       // ... assemble K and F, apply the boundary conditions, solve
       return { c, K, F, quadratureError };   // c: nodal values
     }
   }
   ```

   Contract: `solve(problem, mesh) → { c, K, F, quadratureError }`. The method selector, the charts, the reference
   solution and the convergence study use the new entry automatically.
2. Add convergence tests with known analytical solutions to `tests.js` (`node tests.js`).
3. The solution sheet (`fem_sheet.py`) currently writes only the Galerkin formulation and ignores the `method` field.
   A sheet for a new method needs the symbolic counterpart of the method in the engine; without it the cross-check
   difference on the sheet becomes large and shows the mismatch openly.

## 6. Design rules

- The two kernels stay independent.
- One-way dependencies: each file starts with a section map; a section uses only the sections above it.
- Single state object in the interface: every control changes only `state` and calls `render()`.
- `build_sheet` is the single source of the sheet: the on-screen sheet, the copied text and the `.tex` file come from
  the same blocks.
- Security: `fem_sheet.parse` checks the input token by token against a whitelist; the server listens only on
  `127.0.0.1` and does not accept LaTeX from the client.
- No framework, no build step, no AI calls: every result is computed.
