# Verification, diagnostics and testing

This document describes how FEM Lab checks its own results: well-posedness diagnostics, analytical-solution
validation, the reference solutions used for error measurement, the independent verification layer, the error
measures and their interpretation, and the automated tests. The mathematics is in [MATEMATIK.md](MATEMATIK.md); the
software structure in [MIMARI.md](MIMARI.md).

All numerical thresholds below are **heuristic**. They are relative, have been tested for the cases listed, and are
not universal mathematical constants.

## 1. Two independent kernels

| Kernel | Where | Arithmetic | Used for |
|---|---|---|---|
| Numerical | `index.html` (JavaScript) | 64-bit floating point, composite Gauss quadrature (8 panels × 12 points = 96 points per element) | live charts and summary tiles |
| Symbolic | `fem_sheet.py` (SymPy) | exact / symbolic integrals; numerical only where no closed form exists | solution sheet |

The two kernels are deliberately not connected. Every solution sheet compares them (the **cross-check** shown above
the sheet):

- nodal coefficients $\mathbf{a}$;
- the assembled matrix $\mathbf{K}$ (up to 65 unknowns);
- satisfaction of the prescribed $u$ values;
- the coefficient values parsed by the browser against those parsed by the server (9 points; see §6.1);
- for first-order problems, the residual of the dropped row.

Tolerance: $10^{-9}\times$ scale. If the browser's quadrature-error estimate signals non-smooth coefficients, the
tolerance is relaxed to $10^{-5}$ and the relaxation is stated. The nodal error with respect to the analytical
solution is printed separately; it is the method's error, not a disagreement between kernels.

## 2. Reference solution: analytical or fine mesh

Measuring an error requires a reference. FEM Lab uses one of two, and labels which one:

| Reference | When | What the error means |
|---|---|---|
| **Analytical** | The server finds a closed-form solution with SymPy (`/api/exact`) and it passes validation (§4) | The true error $E = u - u_h$. Tile: "Referans u: analitik". |
| **Fine mesh** | No validated analytical solution, or the page is used without the server | The same method on $n = 128$, $p = 4$. The measures are differences $D = u_{\mathrm{ref}} - u_h$ (labelled $D_{L^2}$, $D_{\max}$, …), not errors. Tile: "Referans u_ref: ince ağ". |

Consequences of the fine-mesh reference:

- It is computed with the same method, so a systematic error of the method **cannot** show up in the charts. The
  independent checks are then the SymPy cross-check, the analytical table on the sheet and the verification layer.
- The reference has its own error, estimated asymptotically as $E_{L^2}(n{=}64, p{=}4)/31$ for $L^2$ and
  $E_{H^1}(n{=}64, p{=}4)/15$ for $H^1$. In the convergence chart this estimate, with a safety factor of 3, is added to
  the measurement floor (hatched region). If a measured $L^2$ difference is smaller than 10 times the estimate, the
  tile is marked unreliable. The estimate is asymptotic and may be optimistic for non-smooth solutions.
- If the coefficient's zero falls on a node of the reference mesh (for example $x^* = 1/2$ for a sign-changing
  $\alpha$), the reference system may be singular; then no error is reported and the application states which system
  (user mesh or reference mesh) is singular.

An analytical solution is used only if the browser's parser can read SymPy's text (solutions containing `erf` or
`asinh`, for example, cannot be read and fall back to the fine mesh). The parsed expression is compared with sample
values sent by the server.

## 3. Well-posedness: continuous and discrete

The application reports the state of a problem at **two separate levels** and never mixes them.

### 3.1 Continuous problem (server, SymPy)

The general solution $u = u_p + \sum_{i=1}^{m} C_i e_i$ (obtained with `dsolve` without initial conditions; $m$ = order)
is substituted into the $q$ conditions that the problem definition imposes, giving a $q \times m$ system
$\mathbf{M}\mathbf{C} = \mathbf{r}$ (Fredholm alternative):

| Ranks | Conclusion |
|---|---|
| $\operatorname{rank}\mathbf{M} = \operatorname{rank}[\mathbf{M}\,\vert \,\mathbf{r}] = m$ | unique solution |
| $\operatorname{rank}[\mathbf{M}\,\vert \,\mathbf{r}] > \operatorname{rank}\mathbf{M}$ | **no solution** |
| $\operatorname{rank}\mathbf{M} = \operatorname{rank}[\mathbf{M}\,\vert \,\mathbf{r}] < m$ | **not unique** |

- Only a **symbolically proven** "no solution" or "not unique" blocks the solution sheet.
- The symbolic rank decision is compared with a 50-digit numerical SVD rank. If they disagree, the result is
  labelled "numerical (not a proof)" and does not block anything.
- If no closed-form general solution exists (most variable-coefficient problems), the continuous problem is reported
  as "could not be classified"; nothing is rejected.
- **Degenerate end point** ($\alpha(x_b) = 0$): classification is done only under an explicit assumption that is
  printed with the result: solutions are bounded on the closed interval (the $C^0$ space of the FEM); the coefficient of
  the fundamental solution that is unbounded at that end is set to zero. Example: $(x u')' = 1$, $u(0) = u(1) = 0$ gives
  "no solution in the class of bounded solutions".
- **Interior zero** of the leading coefficient (including a sign change): no continuous classification is attempted
  ("could not be classified"), because separate constants would be needed on each side.

### 3.2 Discrete system (each kernel separately)

- A condition-number estimate $\hat\kappa_1$ (Hager–Higham) is computed for every solve.
- On pivot breakdown, or if $\hat\kappa > 10^{14}$, an SVD determines rank, null space and compatibility, using the
  normwise backward error $\eta = \|\mathbf{b} - \mathbf{A}\hat{\mathbf{x}}\| / (\|\mathbf{A}\|\,\|\hat{\mathbf{x}}\| + \|\mathbf{b}\|)$
  (Rigal–Gaches).
- For a singular system, compatibility means that $\mathbf{b}$ lies in the range of $\mathbf{A}$: for every
  $\mathbf{y}$ with $\mathbf{y}^{\mathsf T}\mathbf{A} = 0$, $\mathbf{y}^{\mathsf T}\mathbf{b} = 0$. For first-order
  problems the left-null-space test of the full over-determined system is described in
  [MATEMATIK.md §6](MATEMATIK.md#6-first-order-compatibility-left-null-space).
- A singular **discrete** system is reported as a statement about the discrete system only. It is not presented as
  proof that the differential problem has no solution.

### 3.3 Categories

| Category | When |
|---|---|
| K1 No solution (incompatible) | continuous: $\operatorname{rank}[\mathbf{M}\vert \mathbf{r}] > \operatorname{rank}\mathbf{M}$ · discrete: singular and $\eta > 10^{-10}$ |
| K2 Not unique | continuous: $\operatorname{rank}\mathbf{M} = \operatorname{rank}[\mathbf{M}\vert \mathbf{r}] < m$ · discrete: singular and compatible |
| K3 Nearly singular / ill-conditioned | $\hat\kappa > 10^{10}$ (warning) · discrete system singular although the continuous problem is unique (e.g. an initial value problem on a coarse mesh) |
| K4 Degenerate coefficient | $\alpha$ (for first order $\beta$) vanishes or changes sign in the interval: **warning, not rejection** |
| K5 Numerical solver failure | NaN/∞, unclassifiable breakdown |

Heuristic tolerances (relative; tested for scaling by $10^{-8} \dots 10^{8}$ and sizes $N = 3 \dots 257$, where they
gave the same decisions): rank $\sigma \le 10^{-12}\,\sigma_{\max}$; compatibility $\eta \le 10^{-10}$; conditioning
warning at $10^{10}$; numerically singular at $10^{14}$.

A zero or sign change of $\alpha$ alone is not a lack of solutions. A sign-changing $\alpha$ whose zero coincides with a
node can make the discrete system singular; this is reported only as a property of the discrete system.

## 4. Analytical-solution validation

A candidate analytical solution is used only after it passes all of the following:

1. the equation residual $\mathcal{L}(u) - f \equiv 0$;
2. every prescribed $u$ condition;
3. every prescribed $u'$ condition;
4. every prescribed flux condition $\alpha u' = q$, in the form the user entered it (it is not converted into a
   derivative condition);
5. the natural condition only where the formulation adds it to the problem definition (second-order boundary value
   problem, end point without a condition, $\alpha(x_b) \neq 0$); it is not checked at the free end of an initial value
   problem or for first-order problems; at a degenerate free end, boundedness is checked instead;
6. finiteness on the closed interval.

Evidence levels are reported separately: **symbolically verified** (simplification to zero or an exact limit; a
proof) or **numerically verified** (50 digits, $|{\rm value}| \le 10^{-30}\times$ scale; not a proof). The browser
re-checks the conditions with its own parser (tolerance $10^{-9}\times$ scale). A rejected candidate is never used for
error measurement; the reason is shown and the reference falls back to the fine mesh.

## 5. Error measures and convergence

### 5.1 Summary tiles (browser)

| Name | Definition | Note |
|---|---|---|
| $E_{L^2}$ | $\big(\int (u - u_h)^2\,dx\big)^{1/2}$ | true FEM reconstruction, 96 Gauss points per element |
| $\vert e\vert _{H^1}$ | $\big(\int (u' - u_h')^2\,dx\big)^{1/2}$ | seminorm (error of the derivative) |
| $E_{\max}$ | $\max \vert u - u_h\vert $ | 32 equally spaced points per element, then a golden-section search (40 steps) around the largest sample; a numerical estimate (lower bound), not the exact maximum |
| $E_{\text{node}}$ | $\max_i \vert u(x_i) - a_i\vert $ | in 1D the nodes may show superconvergence; do not confuse with $E_{\max}$ |
| Measured rate | $r = \log(E_{\text{coarse}}/E_{\text{fine}}) / \log(h_{\text{coarse}}/h_{\text{fine}})$ | with the true $h$; points at the measurement floor are excluded |

With a fine-mesh reference these are differences ($D$), see §2.

### 5.2 Solution sheet (verification layer)

With a validated analytical solution the sheet reports the nodal error (exact fractions where possible), the pointwise
error (table of midpoints between the first two nodes of each element, and the maximum over 200 samples per element; a
lower bound), $\|u - u_h\|_{L^2}$ and $|u - u_h|_{H^1}$, each in absolute and relative form (relative errors are
"undefined" when the norm of $u$ is zero). Without a validated analytical solution the section states that errors were
not computed.

### 5.3 Interpreting convergence

- Theoretical expectation (smooth solution, second-order boundary value problem, non-vanishing $\alpha$):
  $\|e\|_{L^2} = O(h^{p+1})$ and $|e|_{H^1} = O(h^{p})$. These are asymptotic expectations; the application measures
  rates, it does not impose them.
- For first-order Galerkin there is no guarantee of optimal rates; the $p + 1$ / $p$ rates of initial value problems
  are observed values in the tests, not proven results.
- A rate is computed from the finest consecutive pair whose errors are both above the measurement floor; otherwise it
  is shown as "—", never as 0. Measurement floor: relative roundoff $10^{-11}$ (grows like $1/h$ for $H^1$), plus the
  reference's own error estimate when a fine-mesh reference is used.
- The **convergence study on the sheet** solves the problem with the independent solver on three meshes
  ($m$, $2m$, $4m$ elements, at most 129 unknowns) against the validated analytical solution and reports the observed
  rates, the expected rates (only for the class where theory applies: second-order boundary value problem with
  $\alpha \neq 0$ → $L^2$: $p + 1$, $H^1$: $p$) and an asymptotic-regime criterion: two consecutive rates within $0.2$
  of each other; agreement with the expected rate within $0.3$. Two decreasing errors are not treated as a proven
  rate, and a rate that does not match is reported as a warning, not as a defect of the method (the meshes may be
  pre-asymptotic or the smoothness assumption may fail).
- Strongly convection-dominated problems (e.g. $-0.05u'' + u' = 1$ on a coarse mesh) oscillate; this is known
  behaviour of standard Galerkin, not a bug.

## 6. Independent verification layer (`fem_verify.py`)

Each solution sheet runs these checks. Every item carries a status — passed (`geçti`), failed (`başarısız`),
warning (`uyarı`) or not performed (`yapılmadı`) — on the sheet and above it. A check that could not be performed is
never presented as performed. The module does not call any function of the FEM engine; `tests.py` verifies this from
the source code.

### 6.1 Problem-definition equivalence (validation layer)

Before the engine runs, every input field is parsed by two independent parsers (the token-by-token translator and
SymPy's own grammar) and the results are compared symbolically or at 9 interior points with 30 digits (relative
tolerance $10^{-12}$). The table *input → internal representation → equivalence* opens every sheet. The browser
compares its own (third) parser at the same 9 points. Any disagreement stops the computation.

### 6.2 Checks

| Check | Method | Pass criterion (heuristic) |
|---|---|---|
| Analytical solution | §4 | all conditions hold |
| Operator consistency (patch test) | For $u_p = 1 + t + \dots + t^p$ (exactly representable in $V_h$), $\mathbf{K}\mathbf{a}$ must equal the user's operator applied by direct differentiation and integrated with 24-point Gauss (plus the boundary flux terms); load vector likewise. Detects dropped terms, sign errors, variable-coefficient and assembly errors. | relative gap $\le 10^{-10}$; widened by a quadrature-error estimate for non-smooth data (then status "warning") |
| Independent re-assembly and solve | `verification/ifem.py`: Vandermonde basis, Gauss–Legendre, 30-digit LU; nodal values compared | relative gap $\le 10^{-9}$ (or $\kappa\cdot 10^{-14}$ if larger); only up to 129 unknowns, otherwise "not performed" |
| Solver residual | $\Vert \mathbf{F} - \mathbf{K}\mathbf{a}\Vert _\infty$, backward error $\eta = \Vert \mathbf{F} - \mathbf{K}\mathbf{a}\Vert  / (\Vert \mathbf{K}\Vert \,\Vert \mathbf{a}\Vert  + \Vert \mathbf{F}\Vert )$, $\kappa_\infty$, all in 30 digits | $\eta \le 10^{-12}$; warning above $\kappa = 10^{10}$, failure above $10^{14}$ |
| First-order dropped row | residual $r_d$ recomputed with the independently assembled $\mathbf{K}$, $\mathbf{F}$ | relative gap $\le 10^{-9}$ |
| Distributional residual (second order) | classical residual $R_e$ and jumps $J_k$ from the independent basis; identity $\int \phi_i R_e + \phi_i J = 0$ on interior rows | $\le 10^{-9}$ (widened by the quadrature estimate) |
| Convergence study | §5.3 | see §5.3 |
| Stability indicators | §6.3 | warnings only |

### 6.3 Stability indicators

These are warnings, not rejections. A large value does not invalidate a result; it says that standard Galerkin may
produce oscillations or phase error at this resolution.

| Indicator | Definition | Warning |
|---|---|---|
| Element Péclet number (second order, $\beta \neq 0$) | $\mathrm{Pe} = \max \vert \beta\vert \,h / (2\vert \alpha\vert )$ | $\mathrm{Pe}/p > 1$ (for P1 with constant coefficients $\mathrm{Pe} \le 1$ is the no-oscillation condition; for $p > 1$, $\mathrm{Pe}/p$ is an approximate indicator) |
| Wave resolution (where $\gamma/\alpha > 0$) | $k = \sqrt{\max(\gamma/\alpha)}$, intervals per wavelength $2\pi / (k\,h/p)$ | below 10: phase (pollution) error expected; below 4: the mesh cannot resolve the wave |
| First-order step (where $\beta \neq 0$) | $\lambda = \max\vert \gamma/\beta\vert $, $\lambda h/p$ | above 1: large step; standard Galerkin gives no stability guarantee for this class |

The same definitions and thresholds are used in the browser's diagnostics box.

## 7. Automated tests

| Command | What it checks |
|---|---|
| `node tests.js` | the browser kernel (167 checks): quadrature, shape functions, solutions against analytical results, rates, diagnostics categories, scaling, flux conditions, stability indicators |
| `python3 tests.py` | the server pipeline (221 checks): hand-derived matrices, SymPy ↔ JavaScript ↔ analytical solution, the sheet's structure and notation, analytical-solution validation, diagnostics categories, layer boundaries, LaTeX compilation of every sheet; it also runs `verification/run_verification.py --quick` and `--compare` |
| `python3 verification/run_verification.py --full` | the browser kernel against the independent solver and analytical solutions for benchmarks T1–T16 and edge cases E1–E4b, $n = 4 \dots 32$ |
| `python3 verification/run_corpus.py` | the 80-case input-only corpus (§7a); needs `pip install -r verification/requirements-verify.txt` (scikit-fem) |

`tests.py` skips its LaTeX compilation checks if `pdflatex` is not installed; the CI workflow installs LaTeX so that
they always run.

Permanent benchmark problems in `tests.py` (expected values derived by hand, not taken from the application):

- **A:** $u' + u = 1$, $u(0) = 0$ → $u = 1 - e^{-x}$.
- **B:** $u'' = -2$, $u(0) = u(1) = 0$ → $u = x(1 - x)$ (P1 nodes exact, interior error non-zero, jumps $-2h$).
- **C:** $((1+x)u')' = 4x + 2$, $u(0) = 0$, $(1+x)u'(1) = 4$ (flux condition) → $u = x^2$.

Diagnostic tests cover problems with no solution, infinitely many solutions, pure Neumann compatibility, a
near-resonant but unique problem, a degenerate end point, a sign-changing $\alpha$, solver failure, an under-resolved
Helmholtz problem and a convection-dominated problem.

**Regression baseline:** `verification/baseline/` stores the nodal values of the browser kernel recorded before the
V0.1 changes; `run_verification.py --compare` compares the current values against it (relative tolerance
$10^{-12}$).

**Continuous integration:** `.github/workflows/ci.yml` runs the four commands above on every pull request and on every
push to `main`, with Python 3.9 and 3.12 and Node 22.

## 7a. Input-only adversarial corpus and scikit-fem reference

`verification/corpus/cases.json` holds 80 problem **inputs only** (no expected values): first-order IVPs (10),
second-order BVPs (15), variable coefficients (10), derivative/flux conditions (10), element order and mesh (10),
trigonometric/exponential sources (8), singular/degenerate/compatibility cases (7) and invalid/adversarial inputs (10).
Each case is `{id, meta, input}`; only `input` (the same fields the page sends to `POST /api/sheet`) reaches the
application.

`python3 verification/run_corpus.py` sends every input through the real server pipeline (`server.compute_sheet`) and
the browser kernel (`node tests.js --solve`), and evaluates the results independently:

- invalid inputs must be rejected cleanly by both paths (no crash, no hang, no numerical result);
- valid inputs are compared with **scikit-fem** (`verification/skfem_reference.py`, version pinned in
  `verification/requirements-verify.txt`), which assembles the documented weak form itself with `ElementLineP1`,
  `ElementLineP2` and `ElementLinePp(p)`. Its interior DOFs for $p \ge 3$ are hierarchical, so results are mapped to
  FEM Lab's equispaced nodal basis with $\Phi_{Jm} = \psi_m(x_J)$: $K_{\text{Lagrange}} = \Phi^{-\mathsf T} K\,\Phi^{-1}$.
  Quadrature: Gauss–Legendre of order $2p + 40$, with order $2p + 8$ as an error estimate;
- an analytical solution found by the application is re-checked against the equation and every condition with a separate
  parser, and the $L^2$ / $H^1$ errors written on the sheet are recomputed independently;
- element matrices and load vectors are compared element by element (worst element and entry reported), then the
  global $K$ and the nodal values;
- singular/incompatible verdicts are compared with the rank/conditioning of the scikit-fem system;
- the **continuous** problem is analysed independently of the application's classifier
  (`verification/continuous_problem.py`): (R) regularity; Method 1, linear existence–uniqueness theorem for a complete
  initial-condition set; Method 2, a symbolically verified fundamental system (L[y] = 0, Wronskian ≠ 0) and the boundary
  matrix $B c = r$; Method 3, Fredholm/Green compatibility $\int f v = [\alpha v u' - \alpha u v']_a^b$ for
  $\beta \equiv 0$ and separated conditions; Method 4, direct substitution of a real, continuous candidate. Statuses:
  EXISTENCE / NONEXISTENCE / UNIQUENESS / NONUNIQUENESS _CONFIRMED, COMPATIBILITY_SATISFIED / _VIOLATED,
  INCONCLUSIVE, NOT_APPLICABLE. A continuous verdict is never drawn from a solver failure, from (non)singularity of
  a FEM matrix or from sampling. A confirmed continuous result that contradicts the application's behaviour is a FAIL.

Results (`verification/corpus_results/`: `results.json`, `skfem_results.json`, `continuous_results.json`,
`baseline_comparison.json`, `report.md`) are kept apart from the corpus; `baseline_v1/` is the run before the corrective
changes and the reference of the baseline comparison. Case statuses: PASS, WARN, LIMITATION (warning on an
intentionally challenging case), DIAGNOSED (a valid input rejected with a diagnosis that is consistent with independent
evidence), REJECTED_OK, FAIL, RUNNER_ERROR (a defect of the runner itself, reported separately). Each valid case also
gets a cross-verification matrix row (browser ↔ server, FEM Lab ↔ ifem, scikit-fem element / global / nodal, analytical
solution, continuous existence / uniqueness / compatibility, L²/H¹ recomputation, boundary condition) with
PASS / FAIL / WARN / UNSUPPORTED / INCONCLUSIVE / NOT_RUN; only PASS counts as verified. scikit-fem
uses the same documented weak form, so agreement does not exclude a common-mode error in that formulation; the
analytical solutions and hand-derived element matrices guard against that.

## 8. Limits of the verification

- The independent solver and the convergence study run only up to 129 unknowns (dense 30-digit LU); beyond that they
  are reported as "not performed". The convergence study needs a validated analytical solution.
- For non-smooth coefficients or loads (e.g. `abs`) the Gauss rules are approximate; the patch test and the
  independent solve then report a warning (a quadrature limit, not an error). The browser's quadrature is also
  approximate there and the "quadrature error" tile warns about it.
- The independent implementation was produced in the same development process as the application (the same AI
  coding assistant). It provides *code independence*, not a mathematical proof.
- $E_{\max}$ and the sheet's pointwise maximum are sampling-based lower bounds; a narrow peak can be missed.
