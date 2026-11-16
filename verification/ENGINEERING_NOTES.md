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
