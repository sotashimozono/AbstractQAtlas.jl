# relations/statistics.jl — quantum statistics: occupation functions and
# squeezed-state moment identities.
#
# The occupation functions are FORMS (maps ε ↦ n dispatched on the
# `ParticleStatistics` tag); the squeezed-state moments are closed-form
# identities of the `Squeezed` state family.  Both are generic — no
# model enters anywhere.

"""
    occupation(stat::ParticleStatistics, ε; β=nothing, T=nothing, μ=0) -> Real

Mean occupation number of a single-particle level at energy `ε`:

- `Fermionic()`: Fermi–Dirac `n(ε) = 1 / (e^{β(ε−μ)} + 1)` — bounded in
  `[0, 1]`, particle–hole symmetric about `n(μ) = 1/2`.
- `Bosonic()`: Bose–Einstein `n(ε) = 1 / (e^{β(ε−μ)} − 1)` — requires
  `ε > μ` (throws otherwise: the mean occupation diverges as `ε → μ⁺`,
  and `ε < μ` is unphysical for ideal bosons).

Both reduce to the classical Boltzmann limit
[`boltzmann_occupation`](@ref) for `β(ε − μ) ≫ 1`, and obey the exact
structural identity `n_B(ε) − n_F(ε) = 2 n_B(2ε)` at `μ = 0`.
"""
function occupation(::Fermionic, ε::Real; β=nothing, T=nothing, μ::Real=0)
    b = _beta(; β=β, T=T)
    return 1 / (exp(b * (ε - μ)) + 1)
end
function occupation(::Bosonic, ε::Real; β=nothing, T=nothing, μ::Real=0)
    b = _beta(; β=β, T=T)
    ε > μ || error(
        "Bose–Einstein occupation requires ε > μ (got ε = $ε, μ = $μ): " *
        "the ideal-Bose mean occupation diverges at ε → μ⁺",
    )
    return 1 / (exp(b * (ε - μ)) - 1)
end
export occupation

"""
    boltzmann_occupation(ε; β=nothing, T=nothing, μ=0) -> Real

Classical (Maxwell–Boltzmann) occupation `n(ε) = e^{−β(ε−μ)}` — the
common `β(ε − μ) ≫ 1` limit of both quantum statistics.
"""
function boltzmann_occupation(ε::Real; β=nothing, T=nothing, μ::Real=0)
    b = _beta(; β=β, T=T)
    return exp(-b * (ε - μ))
end
export boltzmann_occupation

"""
    squeezed_variances(r) -> (x=…, p=…)

Quadrature variances of the single-mode squeezed vacuum
[`Squeezed`](@ref)`(r)` along its principal axes (φ = 0 convention,
`ħ = 1`, vacuum variance `1/2`):

`Var(x) = e^{−2r}/2,  Var(p) = e^{+2r}/2`.

The product `Var(x)·Var(p) = 1/4` saturates the Heisenberg bound for
every `r` — squeezing redistributes, never creates, uncertainty.  A
squeezing angle `φ ≠ 0` rotates the principal axes without changing
these eigen-variances.
"""
squeezed_variances(r::Real) = (x=exp(-2r) / 2, p=exp(2r) / 2)
export squeezed_variances

"""
    squeezed_mean_photons(r) -> Real

Mean photon number of the squeezed vacuum, `⟨n⟩ = sinh²(r)` —
equivalently `(Var(x) + Var(p))/2 − 1/2`, since
`⟨x²⟩ + ⟨p²⟩ = 2⟨n⟩ + 1`.
"""
squeezed_mean_photons(r::Real) = sinh(r)^2
export squeezed_mean_photons
