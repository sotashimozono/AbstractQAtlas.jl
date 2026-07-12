# structure/keldysh.jl — the equilibrium Keldysh distribution function.
#
# The one model-independent "form" the Keldysh fluctuation–dissipation
# theorem needs: the distribution function `h(ω)` that locks the Keldysh
# component to the spectral part, `G^K = h(ω)(G^R − G^A)` (see
# `relations/keldysh.jl`, [`KeldyshFDT`](@ref)).  It is fixed by the exchange
# statistics alone:
#
#     h_F(ω) = tanh(βω/2) = 1 − 2 n_F(ω),
#     h_B(ω) = coth(βω/2) = 1 + 2 n_B(ω),
#
# so it is the same object as the occupation functions
# (`relations/statistics.jl`), one algebraic step away — the tie is checked
# in the tests.  Computed directly from `tanh`/`coth` here so this form has no
# load-order dependence on the statistics submodule.

"""
    keldysh_distribution(stat::ParticleStatistics, ω; β=nothing, T=nothing) -> Real

The equilibrium Keldysh distribution function `h(ω)` for exchange statistics
`stat`:

- `Fermionic()`: `h(ω) = tanh(βω/2)` — bounded in `[-1, 1]`,
- `Bosonic()`:   `h(ω) = coth(βω/2)` — diverges as `ω → 0` (the classical
  `2T/ω` limit).

This is the function multiplying the spectral weight in the
fluctuation–dissipation theorem `G^K = h(ω)(G^R − G^A)` ([`KeldyshFDT`](@ref)).
It equals `1 ∓ 2n(ω)` with `n` the Fermi–Dirac/Bose–Einstein
[`occupation`](@ref) at `μ = 0`, and it is odd in `ω` and satisfies
`h(ω) → sign(ω)` (fermions) / `h(ω) → coth` as `T → 0`.

```julia
keldysh_distribution(Fermionic(), 1.0; β=2.0)   # tanh(1.0) ≈ 0.7616
keldysh_distribution(Bosonic(), 1.0; T=0.5)     # coth(1.0) ≈ 1.3130
```
"""
function keldysh_distribution(::Fermionic, ω::Real; β=nothing, T=nothing)
    b = _beta(; β=β, T=T)
    return tanh(b * ω / 2)
end
function keldysh_distribution(::Bosonic, ω::Real; β=nothing, T=nothing)
    b = _beta(; β=β, T=T)
    return coth(b * ω / 2)
end
export keldysh_distribution
