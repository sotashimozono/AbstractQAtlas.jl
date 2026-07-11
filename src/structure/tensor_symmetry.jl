# structure/tensor_symmetry.jl — symmetries of the response tensors.
#
# A nonlinear susceptibility is *essentially* a higher-order tensor, and
# its higher-order structure carries a known, model-independent symmetry:
# **intrinsic permutation symmetry**.  The order-n response
# `χ⁽ⁿ⁾_{α; β₁…βₙ}(ω₁, …, ωₙ)` is invariant under any simultaneous
# permutation of the (field-index, frequency) pairs
# `{(β₁, ω₁), …, (βₙ, ωₙ)}` — a consequence of the fields being
# interchangeable perturbations (Armstrong, Bloembergen, Ducuing &
# Pershan, Phys. Rev. 127, 1918 (1962)).  This layer makes that symmetry
# explicit at the type level, so a consumer knows which components are
# forced equal without re-deriving it.

"""
    intrinsic_permutation_symmetric(quantity) -> Bool
    intrinsic_permutation_symmetric(::Type{<:AbstractQuantity}) -> Bool

Whether the response tensor is invariant under permutation of its field
indices (paired with their frequencies) — the intrinsic permutation
symmetry of nonlinear response.  `true` for the susceptibilities and the
conductivity, `false` otherwise (the default).
"""
function intrinsic_permutation_symmetric(q::AbstractQuantity)
    return intrinsic_permutation_symmetric(typeof(q))
end
intrinsic_permutation_symmetric(::Type{<:AbstractQuantity}) = false
intrinsic_permutation_symmetric(::Type{<:Susceptibility}) = true
intrinsic_permutation_symmetric(::Type{<:DynamicalSusceptibility}) = true
intrinsic_permutation_symmetric(::Type{<:Conductivity}) = true
export intrinsic_permutation_symmetric

"""
    canonical_component(χ) -> typeof(χ)

The canonical representative of `χ`'s intrinsic-permutation-symmetry
class: the same response tensor with its **field** indices sorted (the
response index — the first — is fixed).  Two components are forced equal
by the symmetry iff they share a `canonical_component`
([`permutation_equivalent`](@ref)):

```julia
canonical_component(Susceptibility(:x, :z, :y)) === Susceptibility(:x, :y, :z)
```
"""
function canonical_component(χ::Susceptibility)
    i = indices(χ)
    return Susceptibility(i[1], sort!(collect(i[2:end]))...)
end
function canonical_component(χ::DynamicalSusceptibility)
    i = indices(χ)
    return DynamicalSusceptibility(i[1], sort!(collect(i[2:end]))...)
end
function canonical_component(χ::Conductivity)
    i = indices(χ)
    return Conductivity(i[1], sort!(collect(i[2:end]))...)
end
export canonical_component

"""
    permutation_equivalent(a, b) -> Bool

Whether response tensors `a` and `b` are forced equal by intrinsic
permutation symmetry — i.e. differ only by a permutation of their field
indices.  `χ⁽²⁾_{x;yz} = χ⁽²⁾_{x;zy}` (with the frequencies permuted
accordingly), so they are equivalent; `χ_{x;yz}` and `χ_{y;xz}` are not
(different response index).
"""
function permutation_equivalent(a::AbstractQuantity, b::AbstractQuantity)
    (intrinsic_permutation_symmetric(a) && typeof(a).name === typeof(b).name) ||
        return false
    return canonical_component(a) === canonical_component(b)
end
export permutation_equivalent
