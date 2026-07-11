# core/representations.jl — the space a quantity is expressed in, and its
# Fourier-conjugate.
#
# Physical quantities live in one of two spatial representations
# (real space `r` ↔ momentum space `q`) and/or two temporal ones
# (time `t` ↔ frequency `ω`), related by Fourier transform.  These tags
# make "which space" a first-class property, so the transform structure
# (a static structure factor IS the spatial FT of a real-space
# correlation; `S(q,ω)` IS the space-time FT of `⟨A(r,t)A(0,0)⟩`) is
# explicit.
#
# The CONTINUUM transform lives here as structure; its discrete
# realization on a finite grid is a DFT — the `AbstractFFTs.jl` interface
# (`fft`/`ifft`, `fftfreq` for the conjugate grid) — and belongs to the
# functional sibling (issue #14) together with the grid conventions of
# issue #19.  This package only records WHICH representation each quantity
# is in and WHAT its conjugate is.

"""
    AbstractRepresentation

The space a quantity is expressed in.  Spatial: [`RealSpace`](@ref) `r`
↔ [`MomentumSpace`](@ref) `q`.  Temporal: [`TimeDomain`](@ref) `t` ↔
[`FrequencyDomain`](@ref) `ω`.  Fourier-conjugate spaces are paired by
[`fourier_conjugate`](@ref).
"""
abstract type AbstractRepresentation end
export AbstractRepresentation

"""
    RealSpace <: AbstractRepresentation

Real-space representation (position `r` / site index) — conjugate to
[`MomentumSpace`](@ref) under the spatial Fourier transform.
"""
struct RealSpace <: AbstractRepresentation end
export RealSpace

"""
    MomentumSpace <: AbstractRepresentation

Momentum / reciprocal-space representation (wavevector `q`) — conjugate
to [`RealSpace`](@ref).
"""
struct MomentumSpace <: AbstractRepresentation end
export MomentumSpace

"""
    TimeDomain <: AbstractRepresentation

Time-domain representation (`t`) — conjugate to
[`FrequencyDomain`](@ref) under the temporal Fourier transform.
"""
struct TimeDomain <: AbstractRepresentation end
export TimeDomain

"""
    FrequencyDomain <: AbstractRepresentation

Frequency-domain representation (`ω`) — conjugate to
[`TimeDomain`](@ref).
"""
struct FrequencyDomain <: AbstractRepresentation end
export FrequencyDomain

"""
    fourier_conjugate(rep::AbstractRepresentation) -> AbstractRepresentation

The Fourier-conjugate representation: `RealSpace ↔ MomentumSpace`,
`TimeDomain ↔ FrequencyDomain`.  An involution
(`fourier_conjugate(fourier_conjugate(r)) == r`).
"""
fourier_conjugate(::RealSpace) = MomentumSpace()
fourier_conjugate(::MomentumSpace) = RealSpace()
fourier_conjugate(::TimeDomain) = FrequencyDomain()
fourier_conjugate(::FrequencyDomain) = TimeDomain()
export fourier_conjugate
