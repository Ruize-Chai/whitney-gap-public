# Whitney harmonic-space gap on degenerating triangulations

Numerical reproduction for the manuscript **“Quantitative Failure of Whitney Harmonic Approximation on Degenerating Triangulations”**.

**Author:** Ruize Chai  
**Affiliation:** Huazhong University of Science and Technology, Wuhan, China  
**Status:** Research manuscript (pre-submission); numerical companion repository.

The manuscript constructs periodic zig-zag triangulations of a warped two-torus on which the discrete Whitney complex has the correct first cohomology, but its harmonic one-form space need not converge in the ambient `L²` gap. This repository contains independent finite-element assembly checks for the central finite-mesh formulas.

## Model

The periodic torus `T² = R²/Z²` has metric

```text
g = a(y)^2 dx^2 + dy^2,
a(y) = 1 / (A + Δ cos(2πy)),   A > |Δ| > 0.
```

The zig-zag grid is parametrized by `h = 1/n`, `ε = 1/M`, with `n >= 3` and even `M >= 4`. The two smooth harmonic directions are `dx` and `(A + Δ cos(2πy)) dy`.

For the **unweighted scalar** best-approximation problem, the exact finite-mesh energy is

```text
E(h,ε) = Δ²/2 - 2 Δ² τ_h² κ_h² / Γ_h,
θ = 2πh,
τ_h = 8 sin²(θ/4) / θ²,
κ_h = sin(θ/2) / h,
Γ_h = κ_h² + (1 - cos(θ/2))² / (4 ε²).
```

At `ε/h² → ρ ∈ (0,∞)`, this tends to

```text
Λ_ρ = Δ² π² / [2 (16ρ² + π²)].
```

For the metric-weighted vertical harmonic representative, write `e_h²` for its squared `L²_g` error. The two-dimensional harmonic-space geometry gives the gap identity

```text
gap² = e_h² / (A + e_h²).
```

## Requirements

- Python 3.10+
- NumPy 1.24+
- SciPy 1.10+
- Matplotlib 3.7+ (figure generation)

Install:

```bash
python -m pip install -r requirements.txt
```

## Reproduce the checks

Run from the repository root:

```bash
python reproduction/check_mesh_topology.py
python reproduction/verify_exact_assembly.py
python reproduction/verify_weighted_hodge.py
python reproduction/make_figures.py
```

The three verification scripts print per-grid diagnostics and terminate with a `PASS` line after their assertions succeed. The figure-generation script writes PDF/PNG files to `figures/`.

| Script | Independent check | Default test set |
| --- | --- | --- |
| `check_mesh_topology.py` | Periodic vertex/edge/face counts, Euler characteristic, nondegenerate and unique simplices; also exhibits the excluded `M=2` case | `(n,M)=(3,2),(3,4),(4,4),(5,6),(8,32)` |
| `verify_exact_assembly.py` | Direct triangle-by-triangle *unweighted* P1 stiffness/load assembly against the exact Fourier formula | `(4,8),(6,18),(8,32),(10,50),(12,72)` |
| `verify_weighted_hodge.py` | Direct quadrature of the warped metric, Ritz energy bounds, and principal-angle formula via a separately integrated Gram calculation | `(4,8),(6,18),(8,32),(10,50)` |
| `make_figures.py` | Rebuilds the three-regime scalar-energy plot and the periodic zig-zag mesh illustration | `n=8,…,256` and illustrative `M=8, n=5` |

Default numerical settings are `A=2`, `Δ=0.8` for the metric-weighted checks, `Δ=0.8` for the unweighted check, and 12-point Gauss–Legendre quadrature per transformed triangle coordinate. Plot PDFs use black-and-white-friendly line styles; no display server is required.

### Representative outputs

Unweighted scalar energy (Fourier prediction = direct P1 assembly to floating-point tolerance):

| n | M | `E(h,ε)` |
| ---: | ---: | ---: |
| 4 | 8 | 0.073679946805 |
| 6 | 18 | 0.056836773524 |
| 8 | 32 | 0.050743889851 |
| 10 | 50 | 0.047888112677 |
| 12 | 72 | 0.046327273953 |

Metric-weighted checks with `A=2, Δ=0.8`:

| n | M | `e_h²` | `gap²` |
| ---: | ---: | ---: | ---: |
| 4 | 8 | 0.077187363167 | 0.037159557455 |
| 6 | 18 | 0.068537319130 | 0.033133228246 |
| 8 | 32 | 0.065400855333 | 0.031664969618 |
| 10 | 50 | 0.063929236112 | 0.030974529065 |

The reference runs passed on Python/NumPy/SciPy with scalar formula errors below `2×10^-15` and weighted principal-angle identity differences below `10^-9` (the scripts assert this tolerance). Small platform-dependent floating-point differences are expected.

## Numerical method and scope

The assembly scripts identify periodic nodal degrees of freedom but use unwrapped triangle coordinates to compute the physical gradients. They directly integrate the four triangle families; they do **not** import an implementation of the paper's Fourier reduction. The constant kernel of the scalar stiffness matrix is handled by fixing one nodal degree of freedom.

These programs verify finite-grid instances and provide reproducible diagnostics. The full nonconvergence assertion, asymptotic scaling, and universal validity of the formulas rely on the mathematical proof in the manuscript; finite numerical tests alone do not establish them.

## Citation

Ruize Chai, *Quantitative Failure of Whitney Harmonic Approximation on Degenerating Triangulations*, manuscript in preparation. A bibliographic identifier will be added after public release.
