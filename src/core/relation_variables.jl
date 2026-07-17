# core/relation_variables.jl — the RelationVariable layer.
#
# A relation's variables are keyed by the TYPE of the physical thing they
# denote, not by a formula-letter `Symbol` (which drifts and collides — see
# docs/design/type-keyed-interface.md).  The identity-bearing kinds are:
#
#   AbstractQuantity   observables      (core/quantities.jl)
#   AbstractField      control fields   (core/fields.jl)
#   AbstractCoordinate evaluation point (ω, q)              — here
#   AbstractExponent   critical exps    (α, β, γ, ν, η, …)  — here
#
# and a variable also carries a `support` — WHERE it is evaluated — so the same
# type at different regions/points is a distinct key.  The type-keyed prototype
# uses only the trivial `Global` support; `Region` / point / pair supports land
# with the entanglement layer, but the key is `(type, support)` from day one so
# adding them never forces a re-key (design note R1).

"""
    AbstractCoordinate

Parent for evaluation-coordinate variables — a frequency `ω`, a momentum `q`:
the *point at which* a quantity is evaluated, not a subject of the identity.
Coordinates usually appear as lightweight supplied slots rather than typed keys
(design note R3); the type exists so they *can* be keyed when it matters.
"""
abstract type AbstractCoordinate end
export AbstractCoordinate

"""
    Frequency <: AbstractCoordinate

The frequency `ω` at which a dynamical quantity is evaluated.
"""
struct Frequency <: AbstractCoordinate end
export Frequency

"""
    Momentum <: AbstractCoordinate

The momentum / wavevector `q` at which a quantity is evaluated.
"""
struct Momentum <: AbstractCoordinate end
export Momentum

"""
    AbstractExponent

Parent for critical-exponent variables (`α`, `β`, `γ`, `δ`, `ν`, `η`, `z`).
Typing exponents separates the critical-exponent `β` from the inverse
temperature [`InverseTemperature`](@ref) — the two `:β`s a symbol key conflates.
Concrete exponents are introduced when the criticality domain migrates.
"""
abstract type AbstractExponent end
export AbstractExponent

"""
    RelationVariable

Union of the identity-bearing kinds a relation variable can key on:
[`AbstractQuantity`](@ref), [`AbstractField`](@ref), [`AbstractCoordinate`](@ref),
[`AbstractExponent`](@ref).  The *type* is the variable's identity; the formula
letter in the relation body is a private local binding.
"""
const RelationVariable = Union{
    AbstractQuantity,AbstractField,AbstractCoordinate,AbstractExponent
}
export RelationVariable

# Family erasure: `Susceptibility{I}` → `Susceptibility`, so an index-parametric
# quantity collapses to one node/key.  Defined here (the variable-identity file)
# and shared by the auto-derived `quantities` (relations/interface.jl) and the
# quantity graph (structure/graph.jl).
_family(::Type{T}) where {T} = Base.typename(T).wrapper
_family(q::AbstractQuantity) = _family(typeof(q))

# ─── Support: WHERE a variable is evaluated ─────────────────────────────

"""
    Support

Where / on what a variable is evaluated.  [`Global`](@ref) — the whole system,
no decoration — is the default and the only support the type-keyed prototype
uses; region / point / pair supports arrive with the entanglement layer.
"""
abstract type Support end
export Support

"""
    Global <: Support

The trivial support: a bulk, whole-system quantity with no region/point
decoration.  The default support of every variable.
"""
struct Global <: Support end
export Global
Base.show(io::IO, ::Global) = print(io, "global")

# ─── VariableKey: (type, support) — the collision-proof identity ────────

"""
    VariableKey(type::Type, support::Support = Global())

The identity of a relation variable: the quantity / field / coordinate /
exponent `type` together with its `support`.  This — never a formula-letter
`Symbol` — is what the type-keyed bag, [`relation_report`](@ref), and the
derivation graph match on, so distinct types (and, later, distinct supports of
one type) can never collide.
"""
struct VariableKey
    type::Type
    support::Support
end
VariableKey(t::Type) = VariableKey(t, Global())

Base.:(==)(a::VariableKey, b::VariableKey) = a.type === b.type && a.support == b.support
Base.hash(k::VariableKey, h::UInt) = hash(k.type, hash(k.support, hash(:VariableKey, h)))
function Base.show(io::IO, k::VariableKey)
    print(io, "VariableKey(", nameof(k.type))
    k.support isa Global || print(io, ", ", k.support)
    return print(io, ")")
end
export VariableKey
