# structure/scaling_dimensions.jl — the RG-eigenvalue origin of the
# critical exponents.
#
# The four scaling laws in `relations/scaling.jl` (Rushbrooke, Widom,
# Fisher, Josephson) are written there as CHECKABLE identities among a
# supplied exponent set.  But they are not four independent physical
# facts: they are ALL consequences of a single structural statement — that
# the singular part of the free-energy density is a generalized
# homogeneous function of the reduced temperature `t` and the ordering
# field `h`,
#
#     f_s(t, h) = b^{-d} · f_s(b^{y_t} t, b^{y_h} h)                (∗)
#
# with just TWO relevant RG eigenvalues `y_t`, `y_h` and the spatial
# dimension `d`.  Every equilibrium exponent is a fixed rational function
# of `(y_t, y_h, d)`:
#
#     ν = 1/y_t,  α = 2 − d/y_t,  β = (d − y_h)/y_t,
#     γ = (2y_h − d)/y_t,  δ = y_h/(d − y_h),  η = d + 2 − 2y_h.
#
# Substituting these into the four laws, each collapses to `0` IDENTICALLY
# in `(y_t, y_h, d)` — so `critical_exponents(ScalingDimensions(...))`
# produces, by construction, a set that passes every scaling relation
# exactly.  This is the structural root the axioms are derived from: DECLARE
# ONCE (two eigenvalues + `d`), DERIVE EVERYTHING (the six exponents, hence
# the four laws).  No exponent VALUE is stored here — the eigenvalues are an
# input (measured, or read off an RG fixed point); the package owns only the
# universal `(y_t, y_h, d) ↦ exponents` map.
#
# Hyperscaling caveat: equation (∗) carries the bare `b^{-d}`, so the derived
# set satisfies the hyperscaling law (Josephson `2 − α = dν`, and the `d`
# in `η`) by construction.  Above the upper critical dimension a dangerous
# irrelevant variable spoils (∗) and the true (mean-field) exponents cease
# to obey Josephson except exactly at `d = d_upper`; there the eigenvalue
# parameterization no longer applies (see [`Josephson`](@ref)).

"""
    ScalingDimensions(y_t, y_h, d)

The renormalization-group data of a continuous transition: the thermal and
magnetic relevant eigenvalues `y_t`, `y_h` (the RG-flow exponents of the
reduced temperature and the ordering field) and the spatial dimension `d`.
These are the two-and-a-bit numbers the whole equilibrium exponent set is a
function of, via the homogeneity of the singular free energy
`f_s(t,h) = b^{-d} f_s(b^{y_t}t, b^{y_h}h)` — see
[`critical_exponents`](@ref).

Arguments are promoted to a common type; pass `Rational`s for an exact set
(`ScalingDimensions(1//1, 15//8, 2)` is 2D Ising) or `Float64` for a
numerical fixed point.

```julia
critical_exponents(ScalingDimensions(1//1, 15//8, 2))
# (α = 0//1, β = 1//8, γ = 7//4, δ = 15//1, ν = 1//1, η = 1//4)   ← 2D Ising, exact
```
"""
struct ScalingDimensions{T<:Real}
    y_t::T
    y_h::T
    d::T
end
function ScalingDimensions(y_t::Real, y_h::Real, d::Real)
    return ScalingDimensions(promote(y_t, y_h, d)...)
end
export ScalingDimensions

"""
    critical_exponents(s::ScalingDimensions) -> NamedTuple

The full equilibrium exponent set `(α, β, γ, δ, ν, η)` DERIVED from the
RG eigenvalues in `s` — nothing hand-entered:

- `ν = 1/y_t`                (correlation length, `ξ ∼ |t|^{-ν}`)
- `α = 2 − d/y_t`            (specific heat, `C ∼ |t|^{-α}`)
- `β = (d − y_h)/y_t`        (order parameter, `M ∼ |t|^{+β}`)
- `γ = (2y_h − d)/y_t`       (susceptibility, `χ ∼ |t|^{-γ}`)
- `δ = y_h/(d − y_h)`        (critical isotherm, `M ∼ h^{1/δ}`)
- `η = d + 2 − 2y_h`         (anomalous dimension, `G(r) ∼ r^{-(d-2+η)}`)

With `Rational` eigenvalues the result is exact, and it satisfies
[`Rushbrooke`](@ref), [`Widom`](@ref), [`Fisher`](@ref) and
[`Josephson`](@ref) with residual `≡ 0` for *any* `s`
([`exponents_consistent`](@ref) is `true` by construction).
"""
function critical_exponents(s::ScalingDimensions)
    yt, yh, d = s.y_t, s.y_h, s.d
    return (
        α=2 - d / yt,
        β=(d - yh) / yt,
        γ=(2yh - d) / yt,
        δ=yh / (d - yh),
        ν=1 / yt,
        η=d + 2 - 2yh,
    )
end
export critical_exponents

"""
    critical_exponent(name::Symbol, s::ScalingDimensions) -> Real

A single derived exponent (`:α`, `:β`, `:γ`, `:δ`, `:ν`, or `:η`) of `s`
— a keyed view of [`critical_exponents`](@ref).
"""
critical_exponent(name::Symbol, s::ScalingDimensions) = critical_exponents(s)[name]
export critical_exponent

"""
    scaling_dimensions(; ν, η, d) -> ScalingDimensions

Invert the exponent map: recover the RG eigenvalues from the two
independent exponents that fix them, `y_t = 1/ν` and `y_h = (d + 2 − η)/2`,
at dimension `d`.  Composing with [`critical_exponents`](@ref) closes the
loop — every other exponent (`α, β, γ, δ`) is then reconstructed from just
`(ν, η, d)`, a direct expression of the two-eigenvalue structure:

```julia
s = scaling_dimensions(ν = 1//1, η = 1//4, d = 2)   # 2D Ising eigenvalues
critical_exponents(s).δ                              # 15//1  (δ from ν, η, d alone)
```
"""
function scaling_dimensions(; ν, η, d)
    y_t = 1 / ν
    y_h = (d + 2 - η) / 2
    return ScalingDimensions(y_t, y_h, d)
end
export scaling_dimensions
