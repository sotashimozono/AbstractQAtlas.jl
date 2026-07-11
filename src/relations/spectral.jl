# relations/spectral.jl — the pointwise identities relating the
# dynamical / spectral quantities.
#
# These are the frequency-resolved analogues of the thermodynamic
# identities: each holds at a single (q, ω) point (or, for the sum rule,
# given a supplied frequency integral), so each is a first-class
# `@relation` checkable against numbers.  The transform / BZ-sum / low-ω-
# limit relations that are NOT single-point live in `structure/spectral.jl`
# as the inter-quantity genealogy.
#
# Domain tag :spectral throughout, so `relation_report` can sweep a bag of
# measured propagator/response values and report every identity that
# applies.

"""
    Dyson <: AbstractRelation

The Dyson equation relating the full and bare propagators through the
self-energy, at fixed `(q, ω)`:

`G^{-1} = G₀^{-1} − Σ`.

Scalar (single-band) form; the matrix generalization
`G^{-1} = G₀^{-1} − Σ` in orbital/band space is a tracked follow-up.
Complex-valued — pass complex `G`, `G0`, `Σ`.
"""
@relation :spectral Dyson(G, G0, Σ) = 1 / G - (1 / G0 - Σ)

"""
    SpectralFromGreens <: AbstractRelation

The spectral representation of the retarded Green's function at `(q, ω)`:

`A = −(1/π) Im G^R`   ⟺   `A + Im(G^R)/π = 0`.

Pass `ImGR = Im G^R(q, ω)` (a real number) and the real spectral weight
`A`.
"""
@relation :spectral SpectralFromGreens(A, ImGR) = A + ImGR / π

"""
    SpectralSumRule <: AbstractRelation

The single-band spectral normalization

`∫ A(q, ω) dω = 1`.

Supplied-integral convention: `spectral_integral` is the caller-computed
frequency integral of the spectral function at fixed `q`.
"""
@relation :spectral SpectralSumRule(spectral_integral) = spectral_integral - 1

"""
    DetailedBalance <: AbstractRelation

The finite-temperature detailed-balance condition on the dynamical
structure factor,

`S(q, −ω) = e^{−βω} S(q, ω)`,

a convention-independent consequence of the fluctuation–dissipation
theorem.  Variables: `S_plus = S(q, ω)`, `S_minus = S(q, −ω)`, `ω`, and
`β` (or `T`).
"""
@relation :spectral DetailedBalance(S_plus, S_minus, ω, β) = S_minus - exp(-β * ω) * S_plus

"""
    NMRExponent <: AbstractRelation

The dynamical scaling relation fixing the NMR spin–lattice relaxation
exponent from the operator scaling dimension at a quantum critical
point,

`θ_NMR = 2 Δ_op − 1`   (with `1/T₁ ∝ T^{θ_NMR}` as `T → 0`).

Exact-arithmetic: `Δ_op = 1//8` (1D TFIM QCP) gives `θ_NMR = −3//4`
exactly.
"""
@relation :spectral NMRExponent(θ_NMR, Δ_op) = θ_NMR - (2 * Δ_op - 1)
