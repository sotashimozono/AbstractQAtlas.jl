# structure/spectral.jl — the dynamical-quantity inter-relationship graph.
#
# The frequency-resolved quantities form a chain of *operations*:
#
#   SelfEnergy  --dyson-->  RetardedGreensFunction
#               --neg_im_over_pi-->  SpectralFunction
#               --bz_average-->  DensityOfStates
#
#   DynamicalCorrelation  --spacetime_fourier-->  DynamicalStructureFactor
#   DynamicalSusceptibility  --low_frequency_limit-->  NMRSpinRelaxationRate
#
# Each edge records WHICH operation turns the source quantity into the
# target — the inter-quantity genealogy the response tree
# (structure/response.jl) is for thermodynamics.  Where the operation is
# a single (q, ω)-point identity (`neg_im_over_pi`, `dyson`) the edge
# points at the exact `@relation`; where it is a transform / BZ sum /
# limit (`spacetime_fourier`, `bz_average`, `low_frequency_limit`) there
# is no single-point relation — evaluating it is the functional sibling's
# job (issue #14) — so `origin_relation` returns `nothing`.

"""
    SpectralOrigin(from, via)

One edge of the dynamical-quantity graph: the quantity carrying it is
obtained from quantity type `from` by the operation named `via` — one of
`:dyson`, `:neg_im_over_pi`, `:bz_average`, `:spacetime_fourier`,
`:low_frequency_limit`.
"""
struct SpectralOrigin
    from::Type
    via::Symbol
end

"""
    spectral_origin(quantity) -> Union{SpectralOrigin,Nothing}
    spectral_origin(::Type{<:AbstractQuantity}) -> Union{SpectralOrigin,Nothing}

The dynamical-graph edge of `quantity`: the (`from`, `via`) it is
obtained from, or `nothing` for a source quantity
([`SelfEnergy`](@ref), [`DynamicalCorrelation`](@ref),
[`DynamicalSusceptibility`](@ref)) or a quantity outside the graph.

```julia
spectral_origin(DensityOfStates())      # SpectralOrigin(SpectralFunction, :bz_average)
spectral_origin(SpectralFunction())     # SpectralOrigin(RetardedGreensFunction, :neg_im_over_pi)
spectral_origin(RetardedGreensFunction()) # SpectralOrigin(SelfEnergy, :dyson)
```
"""
spectral_origin(q::AbstractQuantity) = spectral_origin(typeof(q))
spectral_origin(::Type{<:AbstractQuantity}) = nothing

spectral_origin(::Type{RetardedGreensFunction}) = SpectralOrigin(SelfEnergy, :dyson)
function spectral_origin(::Type{SpectralFunction})
    return SpectralOrigin(RetardedGreensFunction, :neg_im_over_pi)
end
spectral_origin(::Type{DensityOfStates}) = SpectralOrigin(SpectralFunction, :bz_average)
function spectral_origin(::Type{DynamicalStructureFactor})
    return SpectralOrigin(DynamicalCorrelation, :spacetime_fourier)
end
function spectral_origin(::Type{NMRSpinRelaxationRate})
    return SpectralOrigin(DynamicalSusceptibility, :low_frequency_limit)
end
export spectral_origin, SpectralOrigin

"""
    origin_relation(via::Symbol) -> Union{AbstractRelation,Nothing}

The exact single-`(q, ω)`-point [`@relation`](@ref) that realizes the
operation `via`, or `nothing` when the operation is a transform / sum /
limit with no pointwise form (its evaluation belongs to the functional
sibling, issue #14):

- `:dyson` → [`Dyson`](@ref),
- `:neg_im_over_pi` → [`SpectralFromGreens`](@ref),
- `:bz_average`, `:spacetime_fourier`, `:low_frequency_limit` → `nothing`.
"""
origin_relation(via::Symbol) = origin_relation(Val(via))
origin_relation(::Val{:dyson}) = Dyson()
origin_relation(::Val{:neg_im_over_pi}) = SpectralFromGreens()
origin_relation(::Val) = nothing
export origin_relation

"""
    spectral_chain(quantity) -> Vector{Any}

The dynamical-graph path from `quantity` back to its source, as the list
of quantity types `[typeof(quantity), from, …, source]`.  A source (or
off-graph) quantity returns the singleton `[typeof(quantity)]`.

```julia
spectral_chain(DensityOfStates())
# [DensityOfStates, SpectralFunction, RetardedGreensFunction, SelfEnergy]
# i.e. ρ ⟵ A ⟵ G^R ⟵ Σ : the density of states is built from the
# self-energy through Dyson, the spectral representation, and the BZ sum.
```
"""
function spectral_chain(q::AbstractQuantity)
    chain = Any[typeof(q)]
    cur = typeof(q)
    for _ in 1:32   # depth cap guards a mis-declared cycle
        o = spectral_origin(cur)
        o === nothing && return chain
        push!(chain, o.from)
        cur = o.from
    end
    return error("spectral_chain: graph exceeded depth cap — cyclic edge?")
end
export spectral_chain
