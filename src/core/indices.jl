# core/indices.jl — internal degrees of freedom (tensor indices).
#
# Many physical quantities are tensors in an internal index space and
# were, in the first cut, blurred into scalars by baking one component
# into the type name (`SusceptibilityZZ`, `MagnetizationX`).  This file
# names the index spaces those tensors live in, and declares the tensor
# traits (`tensor_rank`, `index_spaces`, `indices`) that make a
# quantity's tensor character part of the abstract interface — WITHOUT
# any concrete tensor arithmetic.  A quantity type carries its selected
# index values as type parameters (e.g. `Susceptibility{:x,:y}`), so an
# off-diagonal component is expressible, not just the diagonal ones.

"""
    AbstractIndex

Abstract parent for the internal degrees of freedom a tensor quantity's
indices run over.  Concrete singleton tags — [`SpinAxis`](@ref),
[`SpatialDirection`](@ref), [`OrbitalIndex`](@ref) — name *what* an index
is, so [`index_spaces`](@ref) can report a quantity's index structure
without enumerating the (model-dependent) range of values.
"""
abstract type AbstractIndex end
export AbstractIndex

"""
    SpinAxis <: AbstractIndex

A spin / order-parameter component index `α ∈ {x, y, z, …}` — the index
of magnetizations, susceptibilities, spin correlations, structure
factors.
"""
struct SpinAxis <: AbstractIndex end
export SpinAxis

"""
    SpatialDirection <: AbstractIndex

A spatial / Cartesian direction index `μ` — the index of currents and
transport tensors (e.g. the conductivity `σ_μν`).
"""
struct SpatialDirection <: AbstractIndex end
export SpatialDirection

"""
    OrbitalIndex <: AbstractIndex

An orbital / band / sublattice index — the matrix index of single-
particle propagators (`G_ab`, `Σ_ab`, `A_ab`).  Its range is set by the
model, so the interface names the space without enumerating it.
"""
struct OrbitalIndex <: AbstractIndex end
export OrbitalIndex

# ─── Tensor traits ──────────────────────────────────────────────────────

"""
    tensor_rank(quantity) -> Int
    tensor_rank(::Type{<:AbstractQuantity}) -> Int

The tensor rank of a quantity: `0` for a scalar (energy, specific heat,
partition function, density of states, exponents), `1` for a vector
(magnetization `M_α`), `2` for a rank-2 tensor (susceptibility `χ_αβ`,
conductivity `σ_μν`, propagators `G_ab`, structure factors).  Default
`0`; tensorial families override.
"""
tensor_rank(q::AbstractQuantity) = tensor_rank(typeof(q))
tensor_rank(::Type{<:AbstractQuantity}) = 0
export tensor_rank

"""
    index_spaces(quantity) -> Tuple{Vararg{AbstractIndex}}
    index_spaces(::Type{<:AbstractQuantity}) -> Tuple

The internal index spaces of a quantity, one [`AbstractIndex`](@ref) per
tensor slot (so `length(index_spaces(q)) == tensor_rank(q)`).  Empty for
a scalar; `(SpinAxis(), SpinAxis())` for a susceptibility;
`(OrbitalIndex(), OrbitalIndex())` for a Green's function; etc.
"""
index_spaces(q::AbstractQuantity) = index_spaces(typeof(q))
index_spaces(::Type{<:AbstractQuantity}) = ()
export index_spaces

"""
    indices(quantity) -> Tuple{Vararg{Symbol}}
    indices(::Type{<:AbstractQuantity}) -> Tuple

The *selected* index values of a (component of a) tensor quantity, one
symbol per slot — `indices(Susceptibility(:x, :y)) == (:x, :y)`,
`indices(Magnetization(:z)) == (:z,)`, `indices(Energy()) == ()`.  This
replaces the earlier fused `component` label (`:xx`): the honest form is
one entry per index, so a quantity's component pairing is by the whole
tuple.  Length equals [`tensor_rank`](@ref) for a fully-specified
component.
"""
indices(q::AbstractQuantity) = indices(typeof(q))
indices(::Type{<:AbstractQuantity}) = ()
export indices

"""
    response_order(quantity) -> Int
    response_order(::Type{<:AbstractQuantity}) -> Int

For a response function `χ⁽ⁿ⁾ = ∂ⁿ(output)/∂(field)ⁿ`, the order `n` — the
number of conjugate-field derivatives.  Linear response is `1`; the
second-order nonlinear response is `2`; etc.  A response tensor carries
one output index and `n` field indices, so `n = tensor_rank − 1` for the
response families ([`Susceptibility`](@ref), [`Conductivity`](@ref)).
Returns `0` for quantities that are not response functions (the
default).
"""
response_order(q::AbstractQuantity) = response_order(typeof(q))
response_order(::Type{<:AbstractQuantity}) = 0
export response_order
