# Mathematical formulation

This document describes the mathematics that FEM LAB actually implements: the problem class, the boundary
conditions, the Galerkin and weak formulations, the treatment of first-order and initial-value problems, and the
structure of the generated solution sheet. Verification, well-posedness diagnostics and error measurement are
described in [DOGRULAMA.md](DOGRULAMA.md); the software structure in [MIMARI.md](MIMARI.md).

## 1. Problem class

FEM LAB solves one-dimensional, linear, first- or second-order ordinary differential equations of the form

```math
\mathcal{L}(u) \equiv \big(\alpha(x)\,u'\big)' + \beta(x)\,u' + \gamma(x)\,u = f(x), \qquad x \in \Omega = (x_0, x_1).
```

- $\alpha, \beta, \gamma, f$ are user-entered expressions (polynomial, trigonometric, exponential, …). An empty field
  means $0$.
- If $\alpha \equiv 0$ the equation is first order; $\alpha$ and $\beta$ may not both vanish.
- Accepted syntax: `x`, `pi`, `e`, `+ - * / ^`, `sin cos tan asin acos atan sinh cosh tanh exp ln sqrt abs`, and
  implicit multiplication such as `2x` or `(x+1)(x-1)`. Scientific notation (`1e-8`) is rejected, because `e` is Euler's
  number and `1e-8` would otherwise be read as $1\cdot e - 8$; write `10^-8` instead.

Examples of entering equations:

| Equation | α(x) | β(x) | γ(x) | f(x) | Conditions |
|---|---|---|---|---|---|
| $u'' + u = x$ | `1` | | `1` | `x` | e.g. $u(0) = u(1) = 0$ |
| $-u'' = \sin(\pi x)$ | `-1` | | | `sin(pi*x)` | e.g. $u(0) = u(1) = 0$ |
| $u' + 2u = 0$ | | `1` | `2` | | one condition $u(x_b) = g$ |
| $dy/dx = x^2$, $y(0) = 0$ | | `1` | | `x^2` | `u =` `0` at $x_0$ |
| spring $u'' + u = 0$, $u(0) = 0$, $u'(0) = 1$ on $[0, 10]$ (→ $u = \sin x$) | `1` | | `1` | | `u and u′` at $x_0$: $u = 0$, $u' = 1$; $x_1$ `not given` |

Not supported: nonlinear equations, equations of order higher than two, systems of equations, time-dependent
(partial) differential equations, and more than one space dimension.

## 2. Boundary conditions

Each end point $x_b \in \{x_0, x_1\}$ takes one of the following conditions. The flux at an end point is
$Q = \alpha(x_b)\,u'(x_b)$; $Q_0$ and $Q_1$ denote the fluxes at $x_0$ and $x_1$.

| Type (UI label) | Condition | Role in the formulation |
|---|---|---|
| `u =` | $u(x_b) = g$ | essential; the nodal value is known, the flux $Q$ is unknown and computed after the solve |
| `u′ =` | $u'(x_b) = g$ (derivative condition) | natural; enters the load vector as $Q = \alpha(x_b)\,g$ |
| `α u′ =` | $\alpha(x_b)\,u'(x_b) = q$ (flux condition) | natural; enters the load vector as $Q = q$ |
| `u and u′` | $u(x_b) = g$ and $u'(x_b) = s$ at the same end | initial (or terminal) value problem; see §5 |
| `not given` | none | boundary value problem: natural condition $\alpha u' = 0$; initial value problem: free end |

**Derivative and flux conditions are different and are never converted into each other.** For $\alpha = 1 + x$,
the derivative condition $u'(1) = 4$ means $Q_1 = \alpha(1)\cdot 4 = 8$, while the flux condition $(1+x)\,u'(1) = 4$
means $Q_1 = 4$.

Admissible combinations:

- **Second order, boundary value problem:** the two conditions are given at different end points.
- **Second order, initial/terminal value problem:** `u and u′` at one end, the other end left free. Three conditions
  (for example $u$ and $u'$ at $x_0$ and $u$ at $x_1$) over-determine the problem and are rejected.
- **First order:** exactly one condition of the form $u(x_b) = g$.
- A derivative or flux condition requires $\alpha(x_b) \neq 0$. Otherwise the condition cannot enter the weak form
  (the flux term vanishes) and the input is rejected. At an end point *without* a condition where $\alpha(x_b) = 0$,
  "$\alpha u' = 0$" is not an independent condition (it holds automatically); the sheet states this.
- **Robin conditions** ($a\,u + b\,u' = c$) are **not supported**.

## 3. Discretisation

- $n$ equal elements, $1 \le n \le 64$, element size $h = (x_1 - x_0)/n$.
- Lagrange elements of degree $p = 1, \dots, 4$ (linear, quadratic, cubic, quartic) with equally spaced nodes;
  each element has $p + 1$ nodes and the mesh has $N = np + 1$ nodes, i.e. $N$ nodal coefficients $a_j$. The number of
  **unknowns** is smaller when essential conditions are given: every node with a prescribed value $u(x_b) = g$ is
  removed from the unknowns (e.g. P1, $n = 4$: 5 nodes; 4 unknowns with one Dirichlet end, 3 with two).
- Global numbering: local node $k$ of element $e$ (both 0-based) is global node $e\,p + k$.
- Local coordinate $s = x - x_e \in [0, h]$, so $dx = ds$ (the Jacobian is 1).

Notation used throughout (following Heinrich & Pepper, *The Finite Element Method*: $N$ for shape/interpolation
functions, Ch. 3; $\phi$ for Galerkin weighting functions, Ch. 2):

| Symbol | Meaning |
|---|---|
| $u$, $u_h$ | exact solution (when known), finite element approximation |
| $N_j(x)$ | global shape function of node $j$; on $\Omega_e$ it equals the local shape function $N_k^{(e)}(x)$ with $j = (e-1)p + k$ (1-based), elsewhere 0 |
| $N_k(s)$, $N_k^{(e)}(x)$ | local shape function on the reference element, and on element $e$ |
| $a_j$ | unknown global nodal coefficient; $u_h(x_j) = a_j$ |
| $a^{(e)}_k$ | local coefficient of element $e$; $a^{(e)}_k = a_{(e-1)p+k}$ |
| $c_m$ | coefficient of the element polynomial $u^e(s) = \sum_m c_m s^m$ |
| $R(u_h, x)$ | residual $\mathcal{L}(u_h) - f$ |
| $\phi_i$ | Galerkin weighting function; $\phi_i = N_i$ (Bubnov–Galerkin) |
| $K_{ij}$, $F_i$ | system matrix and load vector entries |

### Shape functions

The shape functions are derived by polynomial interpolation. With $u^e(s) = \sum_{m=0}^{p} c_m s^m$ and nodes
$s_k = (k-1)\,h/p$, the interpolation conditions $u^e(s_k) = a^{(e)}_k$ give the system
$\mathbf{V}\mathbf{c} = \mathbf{a}^{(e)}$ with $V_{km} = s_k^{\,m}$. Solving for $\mathbf{c}$ and regrouping by
$a^{(e)}_k$ gives $u^e(s) = \sum_k a^{(e)}_k N_k(s)$. The construction is checked against an independent
product-formula construction, and the following properties are verified for every degree:

```math
N_i(s_j) = \delta_{ij}, \qquad \sum_k N_k(s) = 1, \qquad \sum_k N_k(s)\, s_k^{\,m} = s^m \quad (m = 0, \dots, p).
```

Neighbouring elements share the coefficient of their common node, so $u_h \in C^0$; $u_h'$ is in general
discontinuous at the nodes.

## 4. Galerkin formulation

With $u_h = \sum_j a_j N_j$ the residual $R(u_h, x) = \mathcal{L}(u_h) - f$ is required to vanish in the weighted
mean

```math
\int_{x_0}^{x_1} \phi_i(x)\, R(u_h, x)\, dx = 0,
```

with weighting functions $\phi_i = N_i$ (Bubnov–Galerkin). Weighting functions associated with a node whose value is
prescribed are removed from the weighting space ($\phi(x_b) = 0$ there).

### 4.1 First order: direct weighted residual (no integration by parts)

For $\alpha \equiv 0$ the residual contains only $u_h'$, which is integrable for a $C^0$ approximation, so no
integration by parts is performed and there is no boundary term:

```math
\int_{x_0}^{x_1} \phi_i\,(\beta u_h' + \gamma u_h)\,dx = \int_{x_0}^{x_1} \phi_i f\,dx
\quad\Longrightarrow\quad
K_{ij} = \int \phi_i\,(\beta N_j' + \gamma N_j)\,dx, \qquad F_i = \int \phi_i f\,dx .
```

This is $(-1)$ times the general second-order formula below; the solution is identical, but the sheet uses the
textbook sign and never labels this formulation a "weak form".

### 4.2 Second order: weak formulation (integration by parts)

For a $C^0$ approximation, $u_h'$ jumps at the nodes, so $(\alpha u_h')'$ is not defined there in the classical sense.
Integration by parts moves one derivative onto the weighting function:

```math
\int_{x_0}^{x_1} \phi_i\,(\alpha u_h')'\,dx = \big[\,\phi_i\,\alpha u_h'\,\big]_{x_0}^{x_1} - \int_{x_0}^{x_1} \phi_i'\,\alpha u_h'\,dx ,
```

which gives the weak form

```math
\int_{x_0}^{x_1} \big(\phi_i'\,\alpha u' - \phi_i\,\beta u' - \phi_i\,\gamma u\big)\,dx = \phi_i(x_1)\,Q_1 - \phi_i(x_0)\,Q_0 - \int_{x_0}^{x_1} \phi_i f\,dx ,
\qquad Q = \alpha u'.
```

Substituting $u_h = \sum_j a_j N_j$ and $\phi_i = N_i$:

```math
\mathbf{K}\mathbf{a} = \mathbf{F}, \qquad
K_{ij} = \int \big(\phi_i'\,\alpha\,N_j' - \phi_i\,\beta\,N_j' - \phi_i\,\gamma\,N_j\big)\,dx, \qquad
\mathbf{F} = \mathbf{F}_b + \mathbf{F}_l ,
```

```math
F_{l,i} = -\int \phi_i f\,dx, \qquad \mathbf{F}_b = [\,-Q_0,\ 0,\ \dots,\ 0,\ Q_1\,]^{\mathsf T}.
```

Row $i$ of the system corresponds to the weighting function $\phi_i$; column $j$ multiplies the unknown $a_j$.
The boundary term at each end is handled as follows:

| Condition at $x_b$ | Boundary term |
|---|---|
| $u(x_b) = g$ | $\phi_i(x_b) = 0$ in the weighting space, so the term is absent from the equations used; $Q$ is computed after the solve from the removed row |
| $u'(x_b) = g$ | $Q = \alpha(x_b)\,g$, known, enters $\mathbf{F}_b$ |
| $\alpha u'(x_b) = q$ | $Q = q$, known, enters $\mathbf{F}_b$ |
| not given (boundary value problem) | $Q = 0$ (natural condition) |
| not given, $\alpha(x_b) = 0$ | $Q = 0$ (degenerate end) |
| free end of an initial value problem | unknown; the weighting function that does not vanish there is removed (§5) |

If $\alpha$ is a negative constant, the element matrices are the negatives of the familiar positive-definite
stiffness matrices; $(-\mathbf{K})\,\mathbf{a} = -\mathbf{F}$ has the same solution, and the sheet shows both forms.

### 4.3 Element integrals, assembly, essential conditions

- Element matrices and load vectors are integrated in the local coordinate $s$. In the solution sheet SymPy
  integrates them exactly whenever a closed form exists; otherwise the integral is computed numerically
  (`mpmath.quad`, adaptive tanh–sinh) and marked with "≈". For constant coefficients and equal element sizes the
  element matrix is the same for every element.
- Assembly: $K_{IJ} = \sum_e K^{(e)}_{ij}$ and $F_I = \sum_e F^{(e)}_i$ with $I = (e-1)p + i$, $J = (e-1)p + j$.
- Essential conditions: the row of a prescribed node is removed (its weighting function is not in the weighting space) and its
  column is moved to the right-hand side, $F_{\mathrm{red},i} = F_i - \sum_{j\ \text{known}} K_{ij} a_j$. The reduced
  system $\mathbf{K}_{\mathrm{red}}\,\mathbf{a}_{\text{unknown}} = \mathbf{F}_{\mathrm{red}}$ is solved.

## 5. Initial value problems — what is proven and what is not

Weak form as in §4.2. If $u(x_0) = u_0$ and $u'(x_0) = s_0$ are given at the left end, then $Q_0 = \alpha(x_0)\,s_0$ is
known and $Q_1$ is unknown. Trial space: $u_h(x_0) = u_0$. Weighting space: the $N_i$ with $\phi(x_1) = 0$ (the last
weighting function is removed). This is a Petrov–Galerkin choice: the weighting space differs from the trial space.

- **Proven:**
  1. Because $\phi_i(x_1) = 0$, the unknown term $Q_1$ drops out.
  2. Trial and weighting spaces have the same dimension $np$, so the system is square.
  3. Consistency: the exact solution satisfies this weak form for every weighting function.
  4. At the continuous level, a solution of this weak form satisfies both the differential equation and
     $u'(x_0) = s_0$ (reverse the integration by parts).
- **At the discrete level:** $u(x_0)$ is satisfied exactly (essential condition); **$u'(x_0)$ only weakly**:
  $u_h'(x_0) \neq s_0$ in general (for example $u'' + u = 0$, $n = 2$, P1 gives $u_h'(0) = 24/25$). The difference
  vanishes as the mesh is refined; the sheet prints both values side by side.
- **Stability is conditional:** for $u'' + \omega^2 u = 0$ and P1 the rows form a three-term recurrence whose roots lie
  on the unit circle only for $\omega h < 2\sqrt{3} \approx 3.46$; above it the discrete solution grows exponentially
  (in the tests an amplitude of $10^9$ at $\omega h = 3.6$). For P2 there is a discrete resonance near
  $\omega^2 h^2 \approx 10$ (the interior node equation vanishes). If the system becomes singular on a coarse mesh,
  the application reports that the mesh is too coarse; the continuous problem itself has a unique solution.
- **Not proven:** stability and convergence in general (variable coefficients, $p \ge 2$). On meshes that resolve
  the solution, an $L^2$ rate of $p + 1$ was measured for five problems with known solutions; this is numerical
  evidence, not a proof. On coarse meshes, phase error accumulates for oscillatory solutions.

## 6. First-order compatibility (left null space)

The first-order weighting space follows the textbook: the weighting function of the node with the prescribed value is removed
(its equation is dropped). The full system, however, $\mathbf{A}\mathbf{a} = \mathbf{b}$ with
$\mathbf{A} \in \mathbb{R}^{N \times (N-1)}$ (all rows, columns of the unknowns), is in general over-determined.
It is compatible exactly when $\mathbf{b}$ lies in the range of $\mathbf{A}$, i.e.

```math
\mathbf{y}^{\mathsf T}\mathbf{A} = 0 \;\Longrightarrow\; \mathbf{y}^{\mathsf T}\mathbf{b} = 0 .
```

Let $d$ be the dropped row and $r_d = \sum_j K_{dj} a_j - F_d$ its residual. If the remaining square block
$\mathbf{A}_2$ is non-singular, the left null space is one-dimensional,
$\mathbf{y} = \big(1 \text{ at } d,\ -\mathbf{A}_2^{-\mathsf T}\mathbf{a}_d\big)$, and $\mathbf{y}^{\mathsf T}\mathbf{b} = -r_d$.

- $r_d = 0$: the full system is compatible; the choice of dropped row does not change the result
  (example: $u' = 2x$).
- $r_d \neq 0$: dropping a different row means a different weighting space and therefore a different solution
  (example: $u' = 3x^2$ on $[0, 2]$, $n = 2$ gives $r = -1$). Which one is more accurate depends on the problem;
  this is not an accuracy proof.

In exact arithmetic $\mathbf{y}$ is computed from the null space of $\mathbf{A}^{\mathsf T}$; in floating point the
scaled measure $|\mathbf{y}^{\mathsf T}\mathbf{b}| / (\|\mathbf{y}\|\,\|\mathbf{b}\|)$ is reported. The compatibility is
never claimed to hold automatically.

## 7. Classical, distributional and weak residual (second order)

For a $C^0$ approximation the residual has three different meanings, and the sheet keeps them apart:

1. **Element-interior (classical) residual** $R_e(x) = \mathcal{L}(u_h) - f$ on $\Omega_e$, computed with classical
   derivatives inside the element. It is *not* the residual of the global function.
2. **Distributional residual.** Because $u_h'$ jumps at the interior nodes $x_k$,

   ```math
   \mathcal{L}(u_h) - f = \sum_e R_e\,\chi_{\Omega_e} + \sum_k J_k\,\delta(x - x_k), \qquad
   J_k = \alpha(x_k)\big[u_h'(x_k^+) - u_h'(x_k^-)\big].
   ```

3. **Weak residual.** The Galerkin equation of an interior weighting function is the identity
   $\int \phi_i R_e\,dx + \sum_k \phi_i(x_k)\,J_k = 0$ (element-wise integration by parts in reverse).

Example: for $u'' = -2$, $u(0) = u(1) = 0$, P1, the classical residual is $R_e = 2$ in every element and the nodal
jumps are $J_k = -2h$; $\int \phi_i R_e\,dx = 2h$ and the two contributions cancel. P1 reproduces the exact nodal
values here, yet the interior error is non-zero.

At a Neumann or flux end, $\alpha u_h'(x_b)$ of the discrete solution generally differs from the prescribed $Q$: the
natural condition is satisfied only weakly.

## 8. Reconstruction

The approximate solution is reconstructed element by element, $u_h^{(e)}(x) = \sum_k a^{(e)}_k N_k^{(e)}(x)$,
with the element's own shape functions (no smoothing). All error measures use this reconstruction.

## 9. Worked example (reference values)

$u' = 2x$, $u(0) = 0$, $x \in (0, 2)$, $n = 2$, P1 (exact solution $u = x^2$):

- $K^{(e)} = \begin{bmatrix} -\tfrac12 & \tfrac12 \\ -\tfrac12 & \tfrac12 \end{bmatrix}$,
  $\mathbf{a} = [\,0,\ 2/3,\ 4\,]^{\mathsf T}$, $u_h = 2x/3$ on $[0, 1]$ and $10x/3 - 8/3$ on $[1, 2]$
  (continuous at $x = 1$: $u_h(1) = 2/3$);
- nodal error $1/3$ (at $x = 1$), $\|e\|_{L^2} = \sqrt{4/135} \approx 0.172133$, $|e|_{H^1} = \sqrt{8/9} \approx 0.942809$;
- P2 and P3 reproduce $x^2$ to machine precision.

These values are checked in `tests.py`.

## 10. Structure of the solution sheet

The sheet contains only headings and mathematics (no prose). Section numbering is dynamic; the weak-form section
exists only for second-order problems.

1. **Given data:** a table *input as typed → internal representation → equivalence check*, the strong form, the
   operator $\mathcal{L}$ and its order, the boundary conditions (labelled essential / derivative / flux / initial),
   $n$, $h$, $p$, $N$, the method, the continuous classification and the status of the analytical solution.
2. Notation.
3. Mesh (local-to-global table).
4. Shape functions (derivation via $\mathbf{V}$, $\delta_{ij}$, partition of unity, completeness, continuity,
   local → global).
5. Approximate solution $u_h = \sum_j a_j N_j$.
6. Galerkin weighted-residual formulation (for first order also the evaluated $\int \phi_i R\,dx = 0$ equations, which
   are checked against the rows of $\mathbf{K}$; for second order the classical and distributional residual).
7. *(second order only)* Weak formulation (integration by parts), with the justification of every boundary term.
8. Element matrices and load vectors (each integral: substituted product → expanded integrand → primitive → value;
   the integration rule is stated).
9. Assembly (contributions to shared entries).
10. Boundary conditions and the reduced system (which rows and columns are removed, and why).
11. Solution ($\mathbf{a}$, end fluxes, first-order compatibility, initial-value slope).
12. Element-wise reconstruction.
13. Error measures (only with a verified analytical solution; otherwise "not computed") and the convergence table.
14. Verification and numerical diagnostics.
15. Warnings and limitations.

Length rules: for $n \le 8$ every element's integrals are written out; for $n > 8$ the first two and the last element
are written out and the others appear in the assembly. Matrices with up to 9 unknowns are printed in full, larger
ones by their first 4 and last 2 rows and columns. Up to 8 free unknowns, and only when the reduced $K$ is rational, the system is solved exactly (fractions);
otherwise in decimals (banded LU with partial pivoting and $\kappa_1$, arithmetically the same as mpmath's dense LU), and the
sheet states the reason next to $\mathbf{a}$. Per-element detail (global $N$, $u_h$, residual) is shown for the first two elements;
per-row equations up to 9 unknowns. If the symbolic computation exceeds 60 s, the sheet is regenerated with numerical
integrals (integration steps are then omitted).

Equally spaced nodes become ill-conditioned for degrees above 4; extending beyond $p = 4$ would require
Gauss–Lobatto nodes.
