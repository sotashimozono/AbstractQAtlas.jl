# structure/graph.jl — the unified quantity-relationship graph.
#
# AbstractQAtlas grows several *separate* typed edges between physical
# quantities, each historically its own ad-hoc accessor:
#
#   * `derivative_edge`          — the response genealogy (χ = ∂M/∂h, …),
#   * `spectral_origin`          — the dynamical graph (A = −ImG/π, …),
#   * `fourier_conjugate_quantity` — conjugate pairs (S(q) ↔ ⟨SS⟩(r)),
#   * the relation ↔ quantity links (`quantities` / `relations_constraining`
#     from quantity_links.jl) — every universal law a quantity must obey.
#
# They share one shape — a *typed edge between quantity kinds* — but there was
# no single queryable graph.  This file unifies them, mirroring the vocabulary
# of QAtlas's MODEL graph (`relations(model) -> Vector{Relation}` with a
# `kind`/`from`/`to`/`detail` edge record): here `related_quantities(q)`
# returns `Vector{QuantityEdge}`, and the same typed-edge core answers "how are
# A and B related?" (`quantity_path`) and exports the whole network
# (`quantity_graph`, `quantity_graph_jsonl`).  Nodes are quantity FAMILIES (the
# index-erased UnionAll — `Susceptibility`, not `Susceptibility{(:z,:z)}`) so
# the structural graph is finite and index variants collapse to one node; the
# concrete index is still used to RESOLVE an edge (χ⁽²⁾ ⟶ χ⁽¹⁾ ⟶ M) before the
# endpoints are normalized.
#
# This is the stdlib-only generic core (issue #57): QAtlas's constraint-edge
# kernel is tightly coupled to its model-graph stores (REGISTRY/REALIZES/…,
# bibkeys, generated coherence checks) and does NOT lift cleanly; the genuinely
# generic part — typed edges + traversal + query + network export — is small
# and lives here, in the base package both atlases share.

"""
    QuantityEdge(kind, from, to, detail)

One typed edge of the quantity-relationship graph, connecting quantity
family `from` to family `to`.  `kind` is one of

- `:derivative` — `to` is the potential/quantity `from` is a field-derivative
  of (the response genealogy, [`derivative_edge`](@ref)); `detail` names the
  field (`"∂/∂h"`, `"∂/∂T"`, …).
- `:spectral` — `from` is obtained from `to` by the dynamical-graph operation
  in `detail` ([`spectral_origin`](@ref)'s `via`: `dyson`, `neg_im_over_pi`, …).
- `:fourier` — `from` and `to` are Fourier conjugates
  ([`fourier_conjugate_quantity`](@ref)); `detail == "fourier"`.
- `:law` — `from` and `to` are co-constrained by a universal relation (they
  appear together in some [`quantities`](@ref)`(rel)`); `detail` names the
  relation (`"SusceptibilityFDT"`, …).

Mirrors QAtlas's `Relation` model-graph edge so the two atlases share one
graph vocabulary (models ⊕ quantities).
"""
struct QuantityEdge
    kind::Symbol
    from::Type
    to::Type
    detail::String
end

function Base.show(io::IO, e::QuantityEdge)
    return print(
        io, "QuantityEdge(:", e.kind, ", ", e.from, " → ", e.to, ", \"", e.detail, "\")"
    )
end

# The index-erased FAMILY of a quantity type: `Susceptibility{(:z,:z)}` and the
# bare `Susceptibility` both normalize to the `Susceptibility` UnionAll, so the
# structural graph has one node per quantity kind.
_family(::Type{T}) where {T} = Base.typename(T).wrapper
_family(q::AbstractQuantity) = _family(typeof(q))

# A short human label for a field type in a derivative edge.
_field_label(::Type{F}) where {F} = string(nameof(F))

# Instantiate a representative of a (possibly parametric) quantity family, or
# `nothing` if none of the standard arities construct — used to probe the
# per-type edge accessors and to seed graph traversal.
function _rep_quantity(T::Type)
    for args in ((), (:x,), (:x, :y))
        try
            return T(args...)
        catch
        end
    end
    return nothing
end

"""
    related_quantities(q) -> Vector{QuantityEdge}

The graph neighborhood of quantity `q` (a quantity instance or type): every
[`QuantityEdge`](@ref) it participates in across ALL edge kinds —
its response-genealogy parent (`:derivative`), its dynamical origin
(`:spectral`), its Fourier conjugate (`:fourier`), and every quantity it is
co-constrained with by a universal law (`:law`).  Endpoints are quantity
families (the index-erased UnionAll).

This is the quantity-graph analogue of QAtlas's `relations(model)`.

```julia
related_quantities(Susceptibility(:z, :z))
# QuantityEdge(:derivative, Susceptibility → Magnetization, "∂/∂MagneticField")
# QuantityEdge(:law, Susceptibility → Magnetization, "SusceptibilityFDT")
# QuantityEdge(:law, Susceptibility → StaticStructureFactor, "StructureFactorSusceptibility")
# …
```
"""
related_quantities(q::AbstractQuantity) = _related_quantities(q, typeof(q))
function related_quantities(T::Type)
    # A concrete type dispatches the per-type accessors directly; a bare family
    # (UnionAll) would fall through to their `<:AbstractQuantity` defaults and
    # miss its index-parametric edge, so instantiate a representative first.
    isconcretetype(T) && return _related_quantities(T, T)
    r = _rep_quantity(T)
    return r === nothing ? _related_quantities(T, T) : _related_quantities(r, typeof(r))
end

# `probe` is a quantity instance OR type accepted by the per-type accessors;
# `ct` is the concrete type used to resolve index-dependent edges.
function _related_quantities(probe, ct::Type)
    edges = QuantityEdge[]
    self = _family(ct)

    e = derivative_edge(probe)
    if e !== nothing
        push!(
            edges,
            QuantityEdge(
                :derivative, self, _family(e.parent), "∂/∂$(_field_label(e.field))"
            ),
        )
    end

    o = spectral_origin(probe)
    if o !== nothing
        push!(edges, QuantityEdge(:spectral, self, _family(o.from), string(o.via)))
    end

    c = fourier_conjugate_quantity(ct)
    if c !== nothing
        push!(edges, QuantityEdge(:fourier, self, _family(c), "fourier"))
    end

    # co-constraint edges: for every relation naming this quantity, link it to
    # each OTHER quantity that relation names (a shared universal law).
    for rel in relations_constraining(probe)
        rname = string(nameof(typeof(rel)))
        for U in quantities(rel)
            fam = _family(U)
            fam === self && continue
            push!(edges, QuantityEdge(:law, self, fam, rname))
        end
    end
    return edges
end

# ── the whole graph: node set + undirected adjacency, built once ──

# recursive concrete leaves of the AbstractQuantity hierarchy
function _quantity_leaves(T::Type=AbstractQuantity, acc::Vector{Type}=Type[])
    for S in subtypes(T)
        isabstracttype(S) ? _quantity_leaves(S, acc) : push!(acc, S)
    end
    return acc
end

const _GRAPH_CACHE = Ref{Union{Nothing,Vector{QuantityEdge}}}(nothing)

"""
    quantity_graph() -> Vector{QuantityEdge}

The entire quantity-relationship graph as a deduplicated list of
[`QuantityEdge`](@ref)s — every `:derivative`, `:spectral`, `:fourier` and
`:law` edge among the constructible quantity families.  Built once from the
type hierarchy and cached.  Edges are undirected in spirit but stored once
per unordered `{from,to}`+`kind`+`detail` (the natural orientation of the
accessor is kept).
"""
function quantity_graph()
    cached = _GRAPH_CACHE[]
    cached === nothing || return cached
    seen = Set{Tuple{Symbol,Type,Type,String}}()
    edges = QuantityEdge[]
    for T in _quantity_leaves()
        q = _rep_quantity(T)
        q === nothing && continue
        for e in related_quantities(q)
            # canonical unordered key so an edge and its mirror collapse
            a, b = e.from, e.to
            key = (e.kind, a, b, e.detail)
            rkey = (e.kind, b, a, e.detail)
            (key in seen || rkey in seen) && continue
            push!(seen, key)
            push!(edges, e)
        end
    end
    _GRAPH_CACHE[] = edges
    return edges
end

"""
    quantity_neighbors(fam) -> Vector{QuantityEdge}

Every edge of [`quantity_graph`](@ref) incident on quantity family `fam`
(as `from` OR `to`), i.e. its neighborhood in the assembled undirected graph
— unlike [`related_quantities`](@ref) this also surfaces edges that name
`fam` as their `to` endpoint (e.g. the quantities whose derivative is `fam`).
"""
function quantity_neighbors(fam::Type)
    f = _family(fam)
    return filter(e -> e.from === f || e.to === f, quantity_graph())
end

"""
    quantity_path(a, b) -> Union{Vector{QuantityEdge},Nothing}

A shortest path in the quantity-relationship graph from quantity `a` to
quantity `b` (instances, families, or types), as the sequence of
[`QuantityEdge`](@ref)s connecting them — the machine answer to "how are `a`
and `b` related?".  `nothing` if they are in different connected components;
the empty vector if `a` and `b` are the same family.

```julia
quantity_path(SpecificHeat(), Magnetization(:z))
# SpecificHeat →∂/∂T Energy →∂/∂β FreeEnergy ←∂/∂h Magnetization   (common root)
```
"""
function quantity_path(a, b)
    src = _family(a isa Type ? a : typeof(a))
    dst = _family(b isa Type ? b : typeof(b))
    src === dst && return QuantityEdge[]
    # undirected adjacency from the assembled graph
    adj = Dict{Type,Vector{QuantityEdge}}()
    for e in quantity_graph()
        push!(get!(adj, e.from, QuantityEdge[]), e)
        push!(get!(adj, e.to, QuantityEdge[]), e)
    end
    # BFS
    prev = Dict{Type,QuantityEdge}()
    visited = Set{Type}((src,))
    frontier = Type[src]
    while !isempty(frontier)
        next = Type[]
        for node in frontier
            for e in get(adj, node, QuantityEdge[])
                other = e.from === node ? e.to : e.from
                other in visited && continue
                push!(visited, other)
                prev[other] = e
                other === dst && return _reconstruct_path(prev, src, dst)
                push!(next, other)
            end
        end
        frontier = next
    end
    return nothing
end

function _reconstruct_path(prev, src::Type, dst::Type)
    path = QuantityEdge[]
    node = dst
    while node !== src
        e = prev[node]
        pushfirst!(path, e)
        node = e.from === node ? e.to : e.from
    end
    return path
end

# ── network export (stdlib-only JSON, mirroring QAtlas's *_jsonl idiom) ──

_json_str(s) = string('"', replace(string(s), '\\' => "\\\\", '"' => "\\\""), '"')

"""
    quantity_graph_jsonl([io=stdout]) -> nothing

Stream the whole [`quantity_graph`](@ref) as JSONL for a network/graph view.
The first line is a summary `{"nodes":N,"edges":M}`; each following line is
one edge object `{"kind":…,"from":…,"to":…,"detail":…}`.  No JSON dependency
— hand-rolled, stdlib-only, matching QAtlas's `relations_jsonl` output shape
so a single consumer can render models ⊕ quantities.
"""
function quantity_graph_jsonl(io::IO=stdout)
    edges = quantity_graph()
    nodes = Set{Type}()
    for e in edges
        push!(nodes, e.from)
        push!(nodes, e.to)
    end
    print(io, "{\"nodes\":", length(nodes), ",\"edges\":", length(edges), "}\n")
    for e in edges
        print(
            io,
            "{\"kind\":",
            _json_str(e.kind),
            ",\"from\":",
            _json_str(nameof(e.from)),
            ",\"to\":",
            _json_str(nameof(e.to)),
            ",\"detail\":",
            _json_str(e.detail),
            "}\n",
        )
    end
    return nothing
end

export QuantityEdge,
    related_quantities,
    quantity_graph,
    quantity_neighbors,
    quantity_path,
    quantity_graph_jsonl
