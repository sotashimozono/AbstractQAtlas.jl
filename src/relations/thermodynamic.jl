# relations/thermodynamic.jl — fluctuation–dissipation identities.
#
# The static fluctuation–dissipation theorem (Callen & Welton, Phys. Rev.
# 83, 34 (1951)) in its thermodynamic form: a response equals β times the
# fluctuation (variance) of its conjugate observable.
#
# One @relation declaration each; the β-or-T keyword convention and all
# solves come from the interface layer (every variable here is affine —
# except β/T themselves, which the generic solve correctly refuses).
#
# Conventions match the QAtlas identities plane (test/identities/):
#
#   SpecificHeat     c_v = β² · Var(E) / N        (E = TOTAL energy)
#   Susceptibility   χ   = β  · Var(M) / N        (M = TOTAL magnetization)
#
# i.e. both responses are PER SITE when `N` is the number of sites and
# the fluctuating quantity is extensive.  `N = 1` (the default) gives
# total responses.  `solve(rel, Val(:C); ...)` / `Val(:χ)` double as the
# estimators downstream Monte-Carlo / ED packages should use, so each
# formula lives in exactly one place.

"""
    SpecificHeatFDT <: AbstractRelation

The energy fluctuation–dissipation identity

`c_v = β² (⟨E²⟩ − ⟨E⟩²) / N = β² Var(E) / N`,

with `E` the total energy and `N` the site count (`N = 1` ⇒ total
specific heat).  Equivalently `c_v = −β² ∂⟨E⟩/∂β / N`: the fluctuation
route and the temperature-response route must agree — that is the
content of the relation, and how it is tested.

```julia
residual(SpecificHeatFDT(); C=c, var_E=v, β=β, N=N)   # c − β²v/N
solve(SpecificHeatFDT(), Val(:C); var_E=v, T=T, N=N)  # the estimator
```
"""
@relation :thermodynamic SpecificHeatFDT(C, var_E, β, N=1) = C - β^2 * var_E / N

"""
    SusceptibilityFDT <: AbstractRelation

The static (zero-frequency) fluctuation–dissipation identity, **per
tensor component**:

`χ_AB = β (⟨M_A M_B⟩ − ⟨M_A⟩⟨M_B⟩) / N = β Cov(M_A, M_B) / N`,

with `M_A` the total (extensive) order-parameter component conjugate to
the field `h_B` and `N` the site count (`N = 1` ⇒ total susceptibility).
Susceptibility is a rank-2 tensor `χ_AB` (`Susceptibility{A,B}`); this
identity relates the `(A, B)` component to the `(A, B)` covariance, so
`var_M` here is `Cov(M_A, M_B)` (the diagonal `A = B` case is the
familiar `Var(M_A)`).  It is the `h → 0` limit of `χ_AB = ∂⟨M_A⟩/∂h_B`;
response route and fluctuation route must agree.
"""
@relation :thermodynamic SusceptibilityFDT(χ, var_M, β, N=1) = χ - β * var_M / N

"""
    LinearResponseFDT <: AbstractRelation

The general static linear-response identity for a perturbation
`H(λ) = H₀ − λ·O` with `[O, H₀] = 0` (classical statistics or a
commuting observable):

`∂⟨O⟩/∂λ = β Var(O)`.

[`SusceptibilityFDT`](@ref) is this relation with `O = M` and a `1/N`
normalization.
"""
@relation :thermodynamic LinearResponseFDT(dO_dλ, var_O, β) = dO_dλ - β * var_O
