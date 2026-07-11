# core/distributions.jl — the "averaged against WHAT" vocabulary.
#
# Statistical distributions / state families as first-class dispatch
# tags carrying their defining parameters, plus particle-statistics
# tags, plus the `ThermalAverage` marker that pairs a quantity with the
# distribution it is averaged in.  This makes the ensemble explicit at
# the type level — instead of an implicit "β kwarg means canonical"
# convention buried in fetch methods.

# ─── Distributions / state families ─────────────────────────────────────

"""
    AbstractDistribution

Abstract parent type for statistical distributions and state families —
the ρ an expectation value is taken against.  Concrete subtypes carry
their defining parameters as typed fields: the three classical ensembles
([`MicroCanonical`](@ref), [`Canonical`](@ref), [`GrandCanonical`](@ref))
and parameterized state families ([`Squeezed`](@ref)).
"""
abstract type AbstractDistribution end
export AbstractDistribution

"""
    MicroCanonical(E; ΔE=0)

Microcanonical ensemble: equal weight on states in the energy window
`|Eᵢ − E| ≤ ΔE/2` (`ΔE = 0` ⇒ exactly-degenerate shell).
"""
struct MicroCanonical{TE<:Real,TW<:Real} <: AbstractDistribution
    E::TE
    ΔE::TW
end
MicroCanonical(E::Real; ΔE::Real=zero(E)) = MicroCanonical(E, ΔE)
export MicroCanonical

"""
    Canonical(β)
    Canonical(; β=nothing, T=nothing)

Canonical (Gibbs) ensemble at inverse temperature `β`, weight
`w(E) = e^{−βE}`.  Constructible from either `β` or `T`.
"""
struct Canonical{TB<:Real} <: AbstractDistribution
    β::TB
end
Canonical(; β=nothing, T=nothing) = Canonical(_beta(; β=β, T=T))
export Canonical

"""
    GrandCanonical(β, μ)
    GrandCanonical(; μ, β=nothing, T=nothing)

Grand-canonical ensemble at inverse temperature `β` and chemical
potential `μ`, weight `w(E, N) = e^{−β(E − μN)}`.
"""
struct GrandCanonical{TB<:Real,TM<:Real} <: AbstractDistribution
    β::TB
    μ::TM
end
GrandCanonical(; μ, β=nothing, T=nothing) = GrandCanonical(_beta(; β=β, T=T), μ)
export GrandCanonical

"""
    Squeezed(r; φ=0.0)

Single-mode squeezed-vacuum state family with squeezing parameter `r`
and squeezing angle `φ` (φ = 0: x-quadrature squeezed).  Not an
ensemble but a parameterized pure-state family — included because its
moments obey generic closed-form identities (see
`relations/statistics.jl`: `squeezed_variances`,
`squeezed_mean_photons`).
"""
struct Squeezed{TR<:Real,TP<:Real} <: AbstractDistribution
    r::TR
    φ::TP
end
Squeezed(r::Real; φ::Real=0.0) = Squeezed(r, φ)
export Squeezed

# ─── Ensemble weights ───────────────────────────────────────────────────

"""
    ensemble_weight(dist, E; N=0) -> Real

The unnormalized statistical weight the distribution assigns to a state
of energy `E` (and particle number `N`, grand-canonical only):

- `MicroCanonical`: indicator of the energy window (0 or 1),
- `Canonical`: `e^{−βE}`,
- `GrandCanonical`: `e^{−β(E − μN)}`.

Normalization is the caller's partition function — summing canonical
weights over a spectrum IS `Z(β)` (cf. [`FreeEnergyFromZ`](@ref)).
"""
ensemble_weight(d::MicroCanonical, E::Real; N=0) = abs(E - d.E) <= d.ΔE / 2 ? 1 : 0
ensemble_weight(d::Canonical, E::Real; N=0) = exp(-d.β * E)
ensemble_weight(d::GrandCanonical, E::Real; N=0) = exp(-d.β * (E - d.μ * N))
export ensemble_weight

# ─── Particle statistics ────────────────────────────────────────────────

"""
    ParticleStatistics

Abstract tag for exchange statistics: [`Fermionic`](@ref) or
[`Bosonic`](@ref).  Used to dispatch statistics-dependent relations
(occupation functions; Wick contractions — the current
[`wick_contraction`](@ref) determinant is the fermionic case, the
bosonic permanent is tracked as a follow-up).
"""
abstract type ParticleStatistics end
export ParticleStatistics

"""
    Fermionic()

Fermi–Dirac exchange statistics tag.
"""
struct Fermionic <: ParticleStatistics end
export Fermionic

"""
    Bosonic()

Bose–Einstein exchange statistics tag.
"""
struct Bosonic <: ParticleStatistics end
export Bosonic

# ─── Thermal-average marker ─────────────────────────────────────────────

"""
    ThermalAverage(quantity, distribution) <: AbstractQuantity

Marker pairing a quantity with the distribution it is averaged in —
`⟨Q⟩_D` at the type level.  Being itself an `AbstractQuantity`, it
composes with the `fetch` verb, so an atlas can make the ensemble an
explicit dispatch axis:

```julia
fetch(model, ThermalAverage(Energy(), Canonical(β)), bc)
```

instead of the implicit "a `β` kwarg means canonical" convention.  The
[`component`](@ref) trait passes through to the wrapped quantity, so
per-component identities keep pairing correctly through the marker.

Fields: `quantity`, `distribution` (accessed directly).
"""
struct ThermalAverage{Q<:AbstractQuantity,D<:AbstractDistribution} <: AbstractQuantity
    quantity::Q
    distribution::D
end
component(::Type{ThermalAverage{Q,D}}) where {Q,D} = component(Q)
export ThermalAverage
