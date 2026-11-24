# Corrective engineering notes — corpus findings (2026-11-16)

Pre-implementation audit (Phase 1), mathematical design (Phase 2) and implementation design (Phase 3) for the defects
found by the 80-case corpus. Evidence after implementation is appended in "Phase 5".
Baseline before any change: `verification/corpus_results/baseline_v1/` (corpus report, scikit-fem results and the raw
server output of all 80 inputs: nodal values, global K, diagnostics, verification items).

## Phase 1 — Execution paths and defects

Execution paths (as traced in the code):

| Path | Call chain |
|---|---|
| Browser → JS kernel | `state` → `compute()` → `buildProblem(problemInputs())` → `makeMesh(x0, x1, n, p)` → `METHODS.galerkin.solve` → `elementSystem` (`REFERENCE[p]`) |
| Server → SymPy engine | `POST /api/sheet` → `server.compute_sheet` → subprocess `fem_sheet.py` → `main` → `make_sheet` → `run_pipeline` → `validate_problem` → `engine_or_diagnosis` (`fem_engine`) → `fem_verify.verify` → `build_sheet` |
| ifem | `fem_verify.verify` → `independent_solution` → `ifem.solve`; `convergence_study` → `ifem.solve` (levels m, 2m, 4m) |
| Analytical solution | `verify` → `classify` → `classify_problem` (thread) → `_classify` (`dsolve`, `num_rank`) → `exact_solution` → `verify_exact` → `error_measures` |
| scikit-fem | `run_corpus.run_case` → `skfem_reference.reference` (own parser, own forms, `asm`, basis transform Φ) |
| Corpus | `run_corpus.main` → `run_case` (server path, `node tests.js --solve`, scikit-fem, independent analytic re-check) → `make_report` |
| CI | `.github/workflows/ci.yml`: tests.py, tests.js, `run_verification.py --full`, `run_corpus.py` |

### D1 — C46: zero free unknowns crash the server (HTTP 500)

- **Input:** `u'' = 1`, n = 1, p = 1, `u(0) = u(1) = 0`: both nodes prescribed, `free = []`, `rows = []`.
- **Engine:** correct (`a = [0, 0]`; `K_ff` is 0×0, `solve_reduced` handles it).
- **Failure points (three, reproduced):**
  1. `fem_verify.independent_solution` → `ifem.solve` → `mp.matrix([])` → `IndexError`.
  2. `fem_verify.convergence_study` → `ifem.solve(prob, 1, 1)` (coarsest level) → the same `IndexError`.
  3. `fem_sheet.build_sheet` → `index_range(tests_used[0] + 1, …)` with an empty weighting space → `IndexError`.
  None is caught; `main` catches only `UserError`, the subprocess dies, `compute_sheet` returns 500.
- **Why tests missed it:** every existing test and benchmark has at least one free unknown.
- **Smallest correction:** (a) `ifem.solve`: no free unknowns → nodal values are the prescribed values, no LU;
  (b) `build_sheet`: an empty weighting space is written explicitly ("no equation to solve"); (c) `fem_verify.verify`:
  every optional check runs in isolation — an exception becomes a *failed* item with the exception text (never PASS),
  and the FEM result and the other checks are still reported.
- **Regression tests:** zero free DOFs (n = 1, both ends `u =`, P1 and P2-free-interior counter-check), sheet produced,
  verification items present; a check that raises is reported as failed; malformed prescribed value is rejected.

### D2 — C76/C77: browser kernel lacks n/p validation

- **Mechanism:** `makeMesh` accepts any n, p. n = 0 → zero elements, `dofs = 1`, both ends fixed → returns `c = [0]`
  (a fabricated "solution"). p = 5 → `REFERENCE[5]` undefined → `TypeError`.
- **Why missed:** UI sliders restrict n ∈ [1, 64], p ∈ [1, 4]; tests.js never passes other values; the server rejects them.
- **Smallest correction:** in the CORE block: `makeMesh` rejects non-integer/≤ 0 n, p ∉ {1,…,4} and a non-finite or
  empty interval (structural soundness; no upper n bound because the internal fine-mesh reference uses n = 128);
  a new `userMesh` additionally enforces n ≤ 64 and is used for every user-supplied mesh (UI, cross-check,
  `tests.js --solve`). Messages identical to the server's (`n ve p tam sayı olmalı.`, `n 1–64, p 1–4 aralığında olmalı.`).
- **Regression tests:** n ∈ {0, −1, 65, 2.5, NaN, missing}, p ∈ {0, 5, 2.5, missing} rejected with `Error` (no
  `TypeError`, no result); n = 1, 64 and p = 1, 4 accepted; internal n = 128 reference still works.

### D3 — C40: complex general solution breaks the numerical rank check

- **Mechanism:** for `((1+x²)u')' = 2` SymPy returns `C2 + (1−C1)log(x−i) + (C1+1)log(x+i)`. Basis
  y₁ = log(x+i) − log(x−i) = 2i(π/2 − atan x) is real up to the factor i, so the boundary matrix M has a purely
  imaginary column (`[[−2i, 0], [iπ/2, 1]]`). Ranks are invariant under column scaling, so the mathematics is fine;
  but `num_rank` calls `mpmath.svd_r` (real only) → `TypeError` → caught → classification "undetermined" → no analytical
  solution, no error norms, no convergence study (silent loss of verification, not a wrong number).
- **Why missed:** no test had a general solution in complex form.
- **Smallest correction:** `num_rank` uses `mpmath.svd_c` for complex matrices (singular values are real; the SVD of A
  uses Aᴴ implicitly), relative threshold unchanged (σ > 10⁻⁴⁰·max(σ_max, 1) at 60 digits). The symbolic rank
  (`Matrix.rank`) is unchanged. The candidate solution is accepted as real only if Im u ≡ 0 is shown symbolically;
  otherwise only after a 50-digit check, labelled numerical (the existing `verify_exact` evidence levels).
- **Regression tests:** real / complex / exactly rank-deficient / nearly rank-deficient matrices; C40 classified
  "unique" with a verified real solution.

### D4 — Element-level scikit-fem comparison missing (verification gap, not an application defect)

- The runner compared global K and nodal values only. **Correction:** the sheet payload carries the engine's element
  matrices and load vectors (float, n ≤ 64 — the ones actually assembled); the runner compares them element by element
  with scikit-fem's single-element assembly mapped to the nodal basis (Φ, see `skfem_reference.py`), reporting the
  worst element and entry.

## Phase 2 — Continuous-level verification: mathematics

Module `verification/continuous_problem.py` (separate from the engine and from the application's own classifier, which
it must not call). Problem: (αu′)′ + βu′ + γu = f on [a, b]; second order iff α ≢ 0.

**Regularity assumptions (R):** α, β, γ, f continuous on the closed interval [a, b] (SymPy `continuous_domain` ⊇ [a, b])
and the leading coefficient (α for second order, β for first order) has no zero on [a, b] (`solveset` on [a, b] returns
the empty set; a `ConditionSet` → INCONCLUSIVE). If (R) cannot be established, Methods 1–3 return INCONCLUSIVE.

**Method 1 — initial value problems** (second order with `u ve u′` at one end; first order with one `u =`):
under (R), linear ODE existence-uniqueness theorem (Picard–Lindelöf for linear systems: global existence and uniqueness
on [a, b]). Conclusions: EXISTENCE_CONFIRMED, UNIQUENESS_CONFIRMED. No computation beyond verifying (R).

**Method 2 — regular two-point BVP via a verified fundamental system:** under (R) the homogeneous equation has a
two-dimensional solution space. A closed form general solution is obtained (SymPy `dsolve`) and *verified independently*:
L[y₁] = L[y₂] = 0 and L[u_p] = f by symbolic simplification, Wronskian W(x*) ≠ 0 at one point (Abel: W never vanishes
if α ≠ 0). Each boundary condition is a linear functional B_k; with c ∈ ℂ²: B c = r, B_kj = B_k[y_j],
r_k = g_k − B_k[u_p]. Conditions: `u =` → u(x_b); `u′ =` → u′(x_b); `α u′ =` → α u′(x_b); `verilmedi` → α u′(x_b) = 0.
- det B ≠ 0 → EXISTENCE_CONFIRMED and UNIQUENESS_CONFIRMED.
- det B = 0 (proved symbolically) → kernel nontrivial: NONUNIQUENESS_CONFIRMED *if* a solution exists;
  rank[B | r] > rank B → NONEXISTENCE_CONFIRMED; = rank B → EXISTENCE_CONFIRMED and NONUNIQUENESS_CONFIRMED.
- det B ≠ 0 decided symbolically, or by a 50-digit evaluation of the exact closed form with |det| ≥ 10⁻³⁰·scale
  (SymPy evalf error control; labelled numerical evaluation of an exact expression). Anything else → INCONCLUSIVE.
- `dsolve` failure, unverifiable pieces, timeouts → INCONCLUSIVE (never nonexistence).

**Method 3 — Fredholm compatibility (self-adjoint case, independent cross-check):** applies only when β ≡ 0
(then L = (αu′)′ + γu is formally self-adjoint) and every condition is separated (`u =`, `u′ =`, `α u′ =`, `verilmedi`),
so the homogeneous boundary operator is self-adjoint. For v in the kernel (L v = 0 with the homogeneous conditions),
Green's identity gives ∫ v L u − ∫ u L v = [α(v u′ − u v′)]ₐᵇ, hence solvability requires

  ∫ₐᵇ f v dx = [α v u′ − α u v′]ₐᵇ  evaluated with the prescribed data,

where at a `u =` end v = 0 (so the unknown u′ drops out) and at the other ends α v′ = 0 (so the unknown u drops out).
Kernel basis from Method 2's verified fundamental system (null vector of B). Equality (symbolic, or 50-digit with the
threshold above) → COMPATIBILITY_SATISFIED; otherwise COMPATIBILITY_VIOLATED → NONEXISTENCE_CONFIRMED.
Not applicable (NOT_APPLICABLE) when β ≢ 0 or a condition is an initial condition.

**Method 4 — direct analytical evidence:** a candidate solution (from Method 2's general solution) substituted into the
ODE and every condition, symbolically. Gives EXISTENCE_CONFIRMED even when (R) fails, never uniqueness.

**Never concluded:** nonexistence from a solver failure; existence or uniqueness from a nonsingular FEM matrix;
nonexistence from a singular FEM matrix; anything from sampling alone. Degenerate/singular problems (α vanishing in
[a, b]) → existence/uniqueness INCONCLUSIVE unless Method 4 produces a verified solution.

## Phase 3 — Implementation design

| Change | Where | Why there |
|---|---|---|
| empty system | `verification/ifem.py` (solve), `fem_sheet.build_sheet` | the failing functions themselves |
| check isolation | `fem_verify.verify` | the only place optional checks are orchestrated |
| complex SVD | `fem_verify._classify.num_rank` | the failing function; symbolic rank untouched |
| real-part acceptance | `fem_verify.exact_solution` | where the candidate is chosen |
| n/p validation | `index.html` CORE (`makeMesh`, new `userMesh`), callers, `tests.js --solve` | CORE is the kernel's input boundary |
| element matrices in payload | `fem_sheet.make_sheet` (payload only) | real pipeline output; engine untouched |
| element comparison, matrix | `verification/run_corpus.py` | runner owns the comparisons |
| continuous verification | new `verification/continuous_problem.py` | must be independent of the application's classifier (`fem_verify`) and of the engine; reuses `skfem_reference.expr` (no third parser) |

The corpus (`verification/corpus/cases.json`) is not changed.

## Phase 5 — Evidence after implementation

Before = `baseline_v1/` (first full run, unchanged engine). After = `corpus_results/` (this run). Corpus unchanged.

| Defect | Reproduction (before) | After | Regression tests |
|---|---|---|---|
| D1 (C46) | `POST /api/sheet` → HTTP 500, `IndexError: list index out of range` (three points, above) | HTTP 200; `a = [0, 0]`; sheet states the empty weighting space; every verification item present, none failed; ifem item *geçti* | `tests.py` "D1: …" (9): zero free DOFs P1; several prescribed values `u(0)=2, u(1)=5`; P2 n = 1 counter-check `a₂ = −1/8` (by hand); malformed values `abc`, `2*(`, `1/0` rejected; a raising check is reported *başarısız* with the exception text, pair → `(item, None)`; `Diagnosis` is not swallowed |
| D2 (C76, C77) | `node tests.js --solve` n = 0 → `c = [0]` (fabricated); p = 5 → `TypeError` | both rejected with Turkish `Error` (same text as the server) | `tests.js` §11 (4): n ∈ {0, −1, 65, 2.5, NaN, missing}, p ∈ {0, 5, 2.5, missing} rejected; n = 1, 64 and p = 1, 4 accepted; internal n = 128 still works; empty/reversed interval rejected |
| D3 (C40) | classification "undetermined" (`svd_r` → `TypeError`), no analytical solution, no error norms | classified unique; real solution `ln(1+x²) + atan x − π/4 − ln 2` verified; L²/H¹ recomputed independently and match the sheet | `tests.py` "D3: …" (8): real full rank, complex full rank (the C40 boundary matrix), exact deficiency real and complex, near deficiency above and below the 10⁻⁴⁰ threshold, complex intermediate with real final rank, C40 end to end |
| D4 | element level not compared | 64 cases compared element by element, worst relative gap 9.2·10⁻¹⁵ | runner check "scikit-fem ↔ FEM Lab eleman matrisleri ve yük vektörleri" |

`numerical_rank` moved from inside `_classify` to module level (above it, one-way dependency kept) so it can be tested;
its arithmetic is unchanged apart from the D3 fix.

**Continuous level** (`continuous_problem.py`, tests `sürekli: …` (9) in `tests.py`): C67 NONEXISTENCE with
C = ∫ f sin πx dx = 1/2 by Method 3, Method 2 consistent; C64 (C = 1) and C66 NONEXISTENCE; C65 EXISTENCE +
NONUNIQUENESS, compatibility satisfied; C46 EXISTENCE + UNIQUENESS; C01 by Method 1; C69/C70 INCONCLUSIVE (no
conclusion without proof).

**Runner defect found during Phase 5 (verification tooling, not the application):** the first post-change run read
C69's discrete message ("… ayrık sistem: … sürekli problem hakkında ayrıca bkz.") as a continuous K1 verdict because
it searched for a substring. Fixed: the category is read from the text before the first `:`; C69 now gets the
"discrete diagnosis ↔ scikit-fem" check again, as in the baseline.

**Baseline comparison:** nodal values of all 63 cases that produced a sheet both times agree to rounding level
(max relative change 6.7·10⁻¹⁶, run-to-run floating-point noise in the decimal solve path; 0 in the first post-change run); the only HTTP change is C46 500 → 200; C76/C77 FAIL → REJECTED_OK; C40 gains the analytical
solution, error norms and convergence study. No new defects, no regressions.

**Observation (not changed):** an empty prescribed value is read as 0 (`parse`: empty text → `'0'`, the same rule as
empty coefficient fields). This is a documented-by-code convention, not a malformed input; changing it would alter UI
behaviour and is outside this task.

## Release-readiness audit (before commit)

Defects found in the new verification code itself (not in the engine) and fixed:

| # | Finding | Fix | Evidence |
|---|---|---|---|
| R1 | Method 2 claimed (Phase 2) a Wronskian check that the code did not perform; a dependent pair would give det B ≡ 0 and a false NONUNIQUENESS/NONEXISTENCE | `fundamental_system` requires W(x*) ≠ 0 at an interior point (symbolic or 50-digit); else INCONCLUSIVE; W is written in the evidence | C68: `W(1/2) ≠ 0 (|W| = 3.13)`; all conclusions unchanged |
| R2 | Method 1 did not check its scope: an over-determined initial set (`u ve u′` + another condition) would have been "EXISTENCE_CONFIRMED" | Method 1 only for exactly one complete initial set; otherwise INCONCLUSIVE | `tests.py` "aşırı belirli koşul takımı …" |
| R3 | det B ≠ 0 without (R) reported compatibility NOT_APPLICABLE "kernel trivial", which is unproven when the fundamental system may not be valid on all of [a, b] | INCONCLUSIVE in that case | C33, C69: NOT_APPLICABLE → INCONCLUSIVE |
| R4 | `real_form` applied the arg → atan identity to Re u only; for C40 Im u ≡ 0 was established only numerically (17 points, 50 digits) | the same exact identity is applied to Im u | C40: Im u ≡ 0 now symbolic (numeric path not entered) |
| R5 | `tests.py` → `continuous_problem` → `skfem_reference` imported scikit-fem, so `tests.py` failed with only `requirements.txt` | shared SymPy reader moved to `verification/ref_parser.py` (still one parser) | `tests.py` imports `continuous_problem` with `skfem` blocked |
| R6 | Phase 5 table said "D1 (8)"; the D1 block runs 9 checks | corrected | `tests.py` 195 → 221: D1 9, D3 8, continuous 9 |

After the fixes: `tests.py` 221/221, `tests.js` 167/167, corpus 80 cases with 0 FAIL, nodal values identical to the
baseline (max relative change 0).

## Phase 6 — CI failures of PR #5 (C48, C50) and silent symbolic-path timeouts (2026-11-16)

CI run 43867191248: Python 3.9 FAIL C50; Python 3.12 FAIL C48 and C50. Both were **HTTP 504** from the server pipeline
(`server.compute_sheet`: symbolic attempt 60 s, then numeric attempt 60 s). No numerical check failed: on 3.9 C48 completed
in 107.5 s and passed every comparison (browser ↔ SymPy 1.1·10⁻¹⁴, scikit-fem elements 1.3·10⁻¹⁵). Failing layer:
**backend execution (performance)**, not formulation, integration, assembly, boundary conditions or solver accuracy. The
3.9/3.12 difference is runner speed around a fixed 60 s limit, not a language-version behaviour.

| # | Finding | Root cause | Fix | Evidence |
|---|---|---|---|---|
| D5 | C48/C50 HTTP 504 | `solve_reduced` (numeric path) built a dense `mpmath.matrix` of the banded K_ff and called `mpmath.cond` (dense inverse) and `mpmath.lu_solve`: O(n³) pure Python. 255 unknowns: 39 s locally (cond 27 s + LU 13 s), i.e. > 60 s on CI runners; the symbolic attempt runs the same solve first | `band_lu` / `band_solve` / `band_lu_solve` / `band_condition` in `fem_sheet.py`: the **same** algorithm as `mpmath.LU_decomp` (scaled partial pivoting, same singularity threshold ‖A‖₁ε, +10 bits for the solve, +20 for the inverse; κ₁ exact via n banded solves, not an estimate), loops restricted to the band (kl, ku + kl after pivoting). Also used for the first-order compatibility solve | 255 unknowns: 39 s → 0.9 s; C50 pipeline 50 s → 13 s. `tests.py`: 40 random banded matrices (incl. zero diagonals → pivoting) agree with `mpmath.lu_solve`/`cond` to < κε, singular → `ZeroDivisionError`; C50 nodes = x(1−x)/2 to 1.2·10⁻¹⁴ in < 20 s |
| D6 | C31 symbolic path > 10 min (server fell back to numeric after 60 s; integral steps lost) | exact (fraction) path chosen by size alone (≤ 8 free unknowns). α = 1/(1+x) puts `log` terms into K; `sp.simplify(det K)` and simplification of the LU solution (≈ 360 ops per entry) never finish | exact path only when the reduced K is **rational** (then det K is an exact rational and the solution is a rational combination of F). Otherwise decimal solve (D5), integrals still exact; reason (`solve_note`) printed next to **a** | C31 > 600 s → 1 s; rational cases unchanged (u'' + u = x: a₂ = −3/44 still exact) |
| D7 | C32 symbolic path 124 s (fallback) | each K entry integrated twice more than needed: the display-only antiderivative used the *expanded* integrand (several × slower for √(1+x)·polynomial), and with β = 0 the integrands N_i'αN_j' and N_j'αN_i' are the same canonical SymPy expression but were integrated separately | antiderivative of the unexpanded integrand (same function; the existing F(h) − F(0) = value check still guards it); per-process memo keyed by (integrand, h, path) — each entry still written out with its own label | C32 124 s → 34 s; test: K^(1)_12 and K^(1)_21 share value and steps, labels distinct |

**Not changed (documented limitations):**
- C35 (57 s) and C63 (50 s) are still close to the 60 s symbolic limit locally (SymPy definite integrals of
  cos(πx)·polynomial and tan x·polynomial). On slower runners they fall back to the numeric sheet (decimal integrals, no
  steps) — a correct but less detailed sheet, not a failure. Further speed-up would need a different integration strategy.
- `classify_problem` / `verify_exact` time out with `thread.join(timeout)`; the abandoned SymPy thread keeps running in the
  same process and competes for the GIL. Results stay correct (a timeout gives INCONCLUSIVE / "undetermined", never a
  verdict), but timing-dependent verdicts differ between runs (C32: "zaman aşımı" on 3.12 vs "not linear in constants" on
  3.9). Killing a SymPy computation needs a subprocess per check; not done here (architecture change).
- p = 5 is not supported by design (equally spaced nodes; Gauss–Lobatto would be required, CLAUDE.md §7).

## Phase 7 — CI failure of PR #7 (C48, C50 again) — verification-layer dense solves (2026-11-24)

CI run 44036363179 (PR #7, branch `hata-grafik-yakinlastirma`, which already contains D5–D7): Python 3.9 FAIL C48 and C50,
both **HTTP 504** after 120.2 s (symbolic 60 s + numeric 60 s); Python 3.12 passed with C48 at 108 s. In PR #6's own run
C48 took 56 s. The same code therefore passes or fails depending on runner load: the cases sit at the limit. Failing layer:
**backend execution (performance)** again, this time in the **verification layer** (Katman 4), which D5 did not touch.

| # | Finding | Root cause | Fix | Evidence |
|---|---|---|---|---|
| D8 | C48 (127 free unknowns) | `fem_verify.solver_check` computed κ∞ with `A ** -1`: dense mpmath inverse, O(n³) pure Python, 7.4 s locally | `ifem.band_inverse`: the **same** arithmetic as `mpmath.inverse` (same LU, +10 bits, unit-vector solves), loops restricted to the band; κ∞ still exact. `SOLVER_LIMIT` unchanged (κ shown up to 129 unknowns, as before) | bit-for-bit equal to `A ** -1` (tests.py) |
| D9 | C48, C50 | `ifem.solve` used dense `mp.lu_solve` (4 calls incl. the convergence study) | `ifem.band_lu_solve`: the same arithmetic as `mpmath.lu_solve` (scaled partial pivoting by row sums, singularity threshold ‖A‖₁ε, +10 bits, same operation order); only operations with exact zeros outside the band are skipped. Implemented in `ifem.py` itself: the verification layer still does not call the engine's `fem_sheet.band_*` | bit-for-bit equal to `mp.lu_solve`, singular → same `ZeroDivisionError` (tests.py, 7 band matrices incl. shifted IVP rows) |
| D10 | C50 (p = 4) | `ifem.basis` re-evaluated N(t), N′(t) (Vandermonde polynomials at 30 digits) for every element and every error sample: 2.9 million `t**m` | values memoised per (t, precision); a copy is returned | — |

Result (local, same machine, `fem_sheet.py` output compared byte for byte with `main` 04c895b): identical for C12, C32,
C35, C40, C46, C48, C49, C50, C54, C62, C63. Numeric path (the one that must finish after a symbolic timeout): C48
10.9 s → 2.8 s, C50 12.0 s → 5.2 s; symbolic path C48 18.3 → 10.5 s, C50 12.5 → 5.4 s, C49 7.0 → 2.4 s.
`tests.py` 228 → 231 (three new checks above). The documented limitations of Phase 6 (C35, C63 near the symbolic limit,
abandoned SymPy threads) are unchanged.
