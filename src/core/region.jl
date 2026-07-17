# core/region.jl — the Region set layer (design doc §5).
#
# A Region is a set of lattice sites — dimension-agnostic (a site is an `Int` in 1D,
# an `NTuple` in ND, or any hashable label), stdlib-only. It is the `support` on
# which a `VonNeumannEntropy` (and, later, a position-dependent field) is evaluated,
# so the entropy inequalities can be auto-discovered over the regions present in a
# bag (`region_report`). The set layer is ND from day one; a geometric layer
# (contiguity, boundary, area law) is a deferred optional lattice extension.

"""
    Region(sites...)
    Region(::AbstractSet)

A subsystem: a set of lattice sites, dimension-agnostic — a site is any hashable
label (`Int` in 1D, `NTuple{D,Int}` in ND, a named block, …).  Supports `∪`, `∩`,
`⊆`, [`disjoint`](@ref), `isempty`, `length`.  The `support` a
[`VonNeumannEntropy`](@ref) is evaluated on; build a region-entropy bag key with
[`entropy`](@ref).

```julia
A, B = Region(1, 2), Region(3, 4)
disjoint(A, B)          # true
A ∪ B                   # Region(1, 2, 3, 4)
Region(1) ⊆ A           # true
```
"""
struct Region{S}
    sites::Set{S}
end
Region(sites...) = Region(Set(sites))
Region(s::AbstractSet) = Region(Set(s))

Base.union(a::Region, b::Region) = Region(union(a.sites, b.sites))
Base.intersect(a::Region, b::Region) = Region(intersect(a.sites, b.sites))
Base.issubset(a::Region, b::Region) = issubset(a.sites, b.sites)
Base.isempty(a::Region) = isempty(a.sites)
Base.length(a::Region) = length(a.sites)
Base.:(==)(a::Region, b::Region) = a.sites == b.sites
Base.hash(a::Region, h::UInt) = hash(a.sites, hash(:Region, h))
function Base.show(io::IO, a::Region)
    return print(io, "Region(", join(sort!(collect(a.sites); by=repr), ", "), ")")
end
export Region

"""
    disjoint(a::Region, b::Region) -> Bool

Whether two regions share no site (`a ∩ b == ∅`) — the precondition of the
bipartite entropy inequalities (subadditivity, Araki–Lieb).
"""
disjoint(a::Region, b::Region) = isempty(intersect(a.sites, b.sites))
export disjoint

"""
    RegionSupport(region::Region) <: Support

The [`Support`](@ref) of a variable evaluated on a [`Region`](@ref):
`VariableKey(VonNeumannEntropy, RegionSupport(A))` keys the entanglement entropy
`S(A)`.  Value-based `==`/`hash` (the `Support` contract), so two content-identical
regions key the same bag entry.
"""
struct RegionSupport{S} <: Support
    region::Region{S}
end
Base.:(==)(a::RegionSupport, b::RegionSupport) = a.region == b.region
Base.hash(a::RegionSupport, h::UInt) = hash(a.region, hash(:RegionSupport, h))
Base.show(io::IO, s::RegionSupport) = print(io, s.region)
export RegionSupport

"""
    entropy(region::Region) -> VariableKey
    entropy(sites...) -> VariableKey

The bag key for the von Neumann entanglement entropy `S(region)` —
`VariableKey(VonNeumannEntropy, RegionSupport(region))`.  Build a region-entropy bag
and auto-discover its inequalities:

```julia
b = bag(entropy(1) => 0.7, entropy(2) => 0.7, entropy(1, 2) => 1.0)   # S(A), S(B), S(A∪B)
region_report(b)                                                       # subadditivity, Araki–Lieb
```
"""
entropy(r::Region) = VariableKey(VonNeumannEntropy, RegionSupport(r))
entropy(sites...) = entropy(Region(sites...))
export entropy
