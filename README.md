# AbstractQAtlas.jl

[![docs: stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://codes.sota-shimozono.com/AbstractQAtlas.jl/stable/)
[![docs: dev](https://img.shields.io/badge/docs-dev-purple.svg)](https://codes.sota-shimozono.com/AbstractQAtlas.jl/dev/)
[![Build Status](https://github.com/QAtlasHub/AbstractQAtlas.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/QAtlasHub/AbstractQAtlas.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Aqua QA](https://raw.githubusercontent.com/JuliaTesting/Aqua.jl/main/badge.svg)](https://github.com/JuliaTesting/Aqua.jl)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

The model-independent core of the QAtlas ecosystem: a registry of physical
quantities and of the universal relations between them.

A quantity's *value* depends on the model; the fact that a susceptibility is a
field-derivative of the magnetization, or that it obeys the
fluctuation–dissipation theorem, does not. AbstractQAtlas owns the second kind
of knowledge — the vocabulary, the domain structure, and the identities that
hold whatever the Hamiltonian — so it can be written down once, tested, and
reused. In the spirit of `AbstractFFTs`, the concrete atlas
([QAtlas.jl](https://github.com/QAtlasHub/QAtlas.jl)) *implements* this package
and supplies the numbers; the dependency never points the other way.

## Quantities, organized by domain

Physical quantities are Julia types, grouped by the domain they belong to:

| domain | quantities |
|---|---|
| **thermodynamics** | `Energy`, `FreeEnergy`, `ThermalEntropy`, `SpecificHeat`, `PartitionFunction`, `Pressure`, … |
| **criticality** | critical exponents `(α, β, γ, δ, ν, η)`, `CorrelationLength`, `MassGap` |
| **correlations** | `RetardedGreensFunction`, `AdvancedGreensFunction`, `KeldyshGreensFunction`, `SpectralFunction`, `DynamicalStructureFactor`, `SelfEnergy` |
| **response / transport** | `Susceptibility{I}` (linear *and* nonlinear), `Conductivity`, `ThermalConductivity`, `Thermopower`, `HallCoefficient` |
| **quantum information** | `VonNeumannEntropy`, `RenyiEntropy`, `MutualInformation`, `Concurrence`, `TopologicalEntanglementEntropy` |
| **quantum foundations** | `KineticEnergy`, `PotentialEnergy`, `EnergyVariance` |
| **topology** | `ChernNumber`, `BerryCurvature`, `BoundaryModeCount` |

Tensor quantities keep their indices — `Susceptibility(:x, :y)` is the
off-diagonal `χ_xy`, `Susceptibility(:x, :y, :z)` the second-order
`χ⁽²⁾_{x;yz}` — rather than collapsing to a scalar.

## Relations between them

The relations are the point. Each is a universal identity among quantities,
written once and carried as a tested, first-class object with three verbs —
`residual` (how far from satisfied), `check` (is it satisfied), and `solve`
(rearrange for one variable):

```julia
using AbstractQAtlas
using AbstractQAtlas: residual, check, solve

# a scaling law, a fluctuation–dissipation identity, a Green's-function equation
solve(Widom(), Val(:δ); β = 1//8, γ = 7//4)              # 15//1, exact
check(SusceptibilityFDT(); χ = 0.5, var_M = 2.0, T = 4.0) # χ = β·Var(M) ?
residual(Dyson(); G = g, G0 = g0, Σ = σ)                 # G⁻¹ = G₀⁻¹ − Σ
```

They span the domains above: the thermodynamic web (`F = U − TS`, the Maxwell
relations, `S = −∂F/∂T`), the critical-exponent scaling laws, the
fluctuation–dissipation theorem in its several forms (thermodynamic, dynamical,
and the Keldysh `G^K = h(ω)(G^R − G^A)`), Kramers–Kronig, Wick's theorem, the
Kubo formula, and the standard topological invariants. Rational inputs give
rational residuals, so an exactly known value satisfies its relation *exactly*,
not merely to floating-point tolerance.

## Following the relations from one quantity to another

Because the relations connect quantities, they form a graph — and that graph is
queryable. You can ask how two quantities are related, or, given what you have
measured, derive another by finding and running a route:

```julia
using AbstractQAtlas: related_quantities, quantity_path, derive

related_quantities(Susceptibility(:z, :z))
#   Magnetization          (as ∂/∂h, and again via fluctuation–dissipation)
#   StaticStructureFactor  (via the structure-factor sum rule)

quantity_path(SpecificHeat(), Magnetization(:z))
#   SpecificHeat — Energy — FreeEnergy — Magnetization   (their common root)

derive(:δ; α = 0//1, β = 1//8, debug = true)
#   DerivationTrace(:δ = 15//1  [indirect])
#     1. Rushbrooke: {α, β} → :γ
#     2. Widom:      {β, γ} → :δ
```

`derive` chains **equalities** only, calls the real `solve` at each step, and —
with `debug = true` — returns the exact route it took, so an indirectly obtained
value is auditable rather than trusted blindly.

## Managed citations

Every relation that rests on the literature carries its source. The references
live in `docs/references.bib`, and CI verifies each DOI against Crossref, so a
citation is either real and resolvable or absent — never fabricated. Kubo (1957)
for the response functions, Callen–Welton (1951) for FDT, Sugiura–Shimizu for
the TPQ construction, and the rest are checked on every commit.

## Installation

```julia
using Pkg
Pkg.add("AbstractQAtlas")
```

Depends only on Julia standard libraries; `ForwardDiff` is a weak dependency
whose loading activates automatic differentiation of the response genealogy.

## Ecosystem

- [QAtlas.jl](https://github.com/QAtlasHub/QAtlas.jl) — the concrete atlas:
  supplies reference values and implements this package's types and relations.

Beyond backing an atlas, AbstractQAtlas is meant as a verification engine: a
quantum-many-body calculation (ED, MPS, TPQ, …) can check its measured
quantities directly against the relation web — fluctuation–dissipation, ensemble
equivalence, scaling, Kubo — instead of re-deriving each identity by hand.
