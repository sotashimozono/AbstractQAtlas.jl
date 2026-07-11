# relations/entanglement.jl — entanglement-entropy relations.
#
# The identities a many-body calculation (MPS/ED) checks its measured
# entanglement against: the Rényi–purity link, the 1D-CFT logarithmic
# growth that reads off the central charge, and Page's average-entropy
# formula for a random pure state.

"""
    RenyiTwoPurity <: AbstractRelation

The Rényi-2 entanglement entropy as (minus) the log purity,

`S_2 = −ln Tr(ρ_A²) = −ln(purity)`

— the `n = 2` member of `S_n = (1−n)⁻¹ ln Tr ρ_A^n`, the one directly
accessible from the [`Purity`](@ref).

Variables: `S2`, `purity`.
"""
@relation :entanglement RenyiTwoPurity(S2, purity) = S2 + log(purity)

"""
    CFTEntanglementSlope <: AbstractRelation

The logarithmic growth of the entanglement entropy of an interval in a
1D conformal field theory reads off the central charge (Calabrese &
Cardy, J. Stat. Mech. (2004) P06002): for a subsystem of size `ℓ`,

`S(ℓ) = (c/3) ln ℓ + const`   ⟹   `dS/d(ln ℓ) = c/3`   (periodic BC;
the open-boundary coefficient is `c/6`).

Supplied-derivative convention: `dS_dlogℓ` is the caller-computed slope
of `S` against `ln ℓ`.  Variables: `dS_dlogℓ`, `c`.
"""
@relation :entanglement CFTEntanglementSlope(dS_dlogℓ, c) = dS_dlogℓ - c / 3

"""
    page_average_entropy(dA, dB) -> Float64

Page's average entanglement entropy of the smaller subsystem `A` for a
Haar-random pure state of a bipartite system `A ⊗ B` with Hilbert-space
dimensions `dA ≤ dB` (Page, Phys. Rev. Lett. 71, 1291 (1993)):

`⟨S_A⟩ = ( Σ_{k=dB+1}^{dA·dB} 1/k ) − (dA − 1)/(2 dB)`.

Nearly maximal, `⟨S_A⟩ ≈ ln dA − dA/(2 dB)`: a random state is almost
maximally entangled, deficit `dA/(2dB)`.  `dA > dB` is symmetric — call
with the smaller dimension first.
"""
function page_average_entropy(dA::Integer, dB::Integer)
    (dA >= 1 && dB >= 1) || error("dimensions must be ≥ 1")
    dA <= dB || return page_average_entropy(dB, dA)   # symmetric; use the smaller
    harmonic = sum(1 / k for k in (dB + 1):(dA * dB))
    return harmonic - (dA - 1) / (2 * dB)
end
export page_average_entropy
