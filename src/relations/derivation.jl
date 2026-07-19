# relations/derivation.jl — the relation registry as a DIRECTED derivation
# graph, and a lazy path-finding solver over it.
#
# Every `@relation` is an exact identity among its variables, and `solve`
# turns it into a computation: given all-but-one variable, produce the last.
# Read that way the whole registry is a DIRECTED graph — a node per relation
# variable (`:f`, `:Z`, `:β`, `:χ`, …), and, for every relation, a directed
# hyperedge `{other variables} →[relation] output` for each variable the
# relation can be solved for.  From a set of KNOWN quantities you can then ask
# what else is reachable, and — lazily — find one route to a target and run it.
#
# Two honesty guarantees, because a chained result is weaker than a directly
# implemented one:
#
#   * `derive` NEVER fabricates an edge: it calls the real `solve`, and a step
#     whose relation is non-affine in the wanted variable simply throws and is
#     skipped (the generic solver refuses non-affine variables by design).  So
#     the reachable set / route reflects what is ACTUALLY computable, not a
#     structural over-estimate.
#   * with `debug=true` the value comes wrapped in a [`DerivationTrace`] naming
#     the exact route (which relation produced each intermediate, from which
#     inputs) and flagged `indirect` — an indirectly derived number carries its
#     provenance so it is never mistaken for a directly-verified one.

"""
    DerivationStep(relation, output, inputs)

One directed edge of the derivation graph: `relation` computes the variable
`output::Symbol` from the variables `inputs::Tuple{Vararg{Symbol}}` (its other
variables), via [`solve`](@ref)`(relation, Val(output); inputs...)`.
"""
struct DerivationStep
    relation::AbstractRelation
    output::Symbol
    inputs::Tuple{Vararg{Symbol}}
end

function Base.show(io::IO, s::DerivationStep)
    ins = join(string.(s.inputs), ", ")
    return print(io, nameof(typeof(s.relation)), ": {", ins, "} → :", s.output)
end

"""
    DerivationTrace

The meta-information returned by [`derive`](@ref)`(...; debug=true)`: the
`target` symbol, its computed `value`, the ordered `steps` that produced it,
and `indirect` — `false` only when the target was among the supplied knowns
(a direct value), `true` when it was derived through one or more relations.
The trace exists for SAFETY: an indirectly derived value is auditable by its
route rather than trusted blindly.
"""
struct DerivationTrace
    target::Symbol
    value::Any
    steps::Vector{DerivationStep}
    indirect::Bool
end

function Base.show(io::IO, t::DerivationTrace)
    tag = t.indirect ? "indirect" : "direct"
    println(io, "DerivationTrace(:", t.target, " = ", t.value, "  [", tag, "])")
    if isempty(t.steps)
        print(io, "  (given directly)")
    else
        for (i, s) in enumerate(t.steps)
            print(io, "  ", i, ". ", s)
            i < length(t.steps) && println(io)
        end
    end
end

const _DERIV_STEPS = Ref{Union{Nothing,Vector{DerivationStep}}}(nothing)

"""
    derivation_steps() -> Vector{DerivationStep}

Every candidate directed edge of the derivation graph: for each registered
relation and each of its variables, the edge that would compute that variable
from the others.  These are STRUCTURAL candidates — a relation appears as an
edge for a variable it is not affine in too; [`derive`](@ref) is the honest
evaluator that discovers, by actually calling [`solve`](@ref), which edges
fire.  Built once and cached.
"""
function derivation_steps()
    cached = _DERIV_STEPS[]
    cached === nothing || return cached
    steps = DerivationStep[]
    for rel in all_relations()
        # EQUALITIES ONLY: `solve` on an inequality returns the SATURATION
        # (the tight bound), not an equational derivation — letting that into
        # the graph would derive a bound and present it as a computed value.
        rel isa AbstractInequality && continue
        vs = variables(rel)
        for out in vs
            ins = Tuple(v for v in vs if v !== out)
            push!(steps, DerivationStep(rel, out, ins))
        end
    end
    _DERIV_STEPS[] = steps
    return steps
end

# Try to fire a step against a value dict; return the solved value or `nothing`
# (the relation is non-affine in `output`, or the solve is otherwise refused).
function _try_step(step::DerivationStep, known::AbstractDict)
    all(v -> haskey(known, v), step.inputs) || return nothing
    try
        return solve(
            step.relation, Val(step.output); (v => known[v] for v in step.inputs)...
        )
    catch
        return nothing
    end
end

# Forward-chaining closure: keep firing any step whose inputs are all known
# until nothing new is produced.  Records the ordered steps actually used.
# Stops early once `stop` (if given) becomes known.
function _forward_chain(known::Dict{Symbol,Any}, stop::Union{Symbol,Nothing})
    used = DerivationStep[]
    progress = true
    while progress && !(stop !== nothing && haskey(known, stop))
        progress = false
        for step in derivation_steps()
            haskey(known, step.output) && continue
            v = _try_step(step, known)
            v === nothing && continue
            known[step.output] = v
            push!(used, step)
            progress = true
            stop !== nothing && haskey(known, stop) && break
        end
    end
    return used
end

# Keep only the steps that actually contribute to `target`, in dependency order.
function _prune(used::Vector{DerivationStep}, target::Symbol, given::Set{Symbol})
    by_output = Dict(s.output => s for s in used)
    needed = DerivationStep[]
    seen = Set{Symbol}()
    function visit(sym)
        (sym in given || sym in seen) && return nothing
        push!(seen, sym)
        step = get(by_output, sym, nothing)
        step === nothing && return nothing
        for inp in step.inputs
            visit(inp)
        end
        return push!(needed, step)   # post-order ⇒ inputs precede their consumer
    end
    visit(target)
    return needed
end

"""
    derivable(; knowns...) -> Set{Symbol}

The set of quantity variables COMPUTABLE from the supplied known values — the
honest reachability of the derivation graph.  Each element is a variable that
can be obtained, directly or through a chain of [`solve`](@ref)s, from the
knowns; the knowns themselves are included.

```julia
derivable(; Z = 2.0, β = 1.0)          # ⊇ Set([:Z, :β, :f])  — F = −β⁻¹ln Z reachable
```

Honest, not structural: a symbol appears only if some relation is actually
affine-solvable for it along the way (non-affine steps are skipped, exactly as
[`derive`](@ref) would skip them).
"""
function derivable(; knowns...)
    known = Dict{Symbol,Any}(pairs(knowns))
    _forward_chain(known, nothing)
    return Set(keys(known))
end

"""
    derive(target::Symbol; debug=false, knowns...) -> value | DerivationTrace

Lazily derive `target` from the supplied known values by finding ONE route
through the derivation graph and running it.  Returns the computed value; with
`debug=true` returns a [`DerivationTrace`](@ref) instead — the value plus the
exact route (which relation produced each intermediate, from which inputs) and
whether it was `indirect`.

```julia
derive(:f; Z = 2.0, β = 1.0)                 # -0.6931…   (F = −β⁻¹ ln Z)
derive(:f; Z = 2.0, β = 1.0, debug = true)   # DerivationTrace: 1. FreeEnergyFromZ: {Z, β} → :f
```

The route is discovered by forward chaining with the REAL [`solve`](@ref), so
a step whose relation is non-affine in its output is skipped, not faked.
Throws if `target` is not reachable from the knowns.
"""
function derive(target::Symbol; debug::Bool=false, knowns...)
    known = Dict{Symbol,Any}(pairs(knowns))
    given = Set(keys(known))
    if haskey(known, target)
        return if debug
            DerivationTrace(target, known[target], DerivationStep[], false)
        else
            known[target]
        end
    end
    used = _forward_chain(known, target)
    haskey(known, target) || error(
        "derive: :$target is not reachable from $(collect(given)) — no chain of " *
        "affine-solvable relations connects them (some steps may need a supplied " *
        "derivative or a value the knowns do not provide).",
    )
    debug || return known[target]
    return DerivationTrace(target, known[target], _prune(used, target, given), true)
end

"""
    derivation_graph() -> KnowledgeGraph{Symbol}

The derivation graph as a [`KnowledgeGraph`](@ref) instance, for a network VIEW
and structural inspection: one DIRECTED edge `input →[relation] output` per
(step, input), the simple-graph projection of the [`derivation_steps`](@ref)
hyperedges.

!!! warning "Structural, not computational"
    Each `@relation`'s output needs ALL of its inputs, but the projection fans
    each hyperedge out to one edge per input — so [`graph_reachable`](@ref) on
    this graph OVER-APPROXIMATES computability (a single known input already
    "reaches" the output, and non-affine outputs are edges too).  For the
    HONEST "can I actually compute it" use [`derivable`](@ref) / [`derive`](@ref),
    which require every input and call the real `solve`.  Use this graph for
    rendering and structural connectivity only.

Edge `kind` is the relation's name.
"""
function derivation_graph()
    edges = TypedEdge{Symbol}[]
    for step in derivation_steps()
        rname = Symbol(nameof(typeof(step.relation)))
        for inp in step.inputs
            push!(
                edges,
                TypedEdge(
                    rname, inp, step.output, string(nameof(typeof(step.relation))), true
                ),
            )
        end
    end
    return KnowledgeGraph(edges)
end

# ─── Type-keyed derivation (collision-proof) ────────────────────────────
#
# The symbol-keyed derivation above keys on private formula LETTERS, which COLLIDE
# across relations (`:S` = thermal entropy / thermopower / structure factor / vN
# entropy): `derive(:Π; S = entropy, T)` silently fires `KelvinRelation`, reusing an
# entropy as a Seebeck coefficient (the demonstrated bug the whole type-keyed
# redesign exists to kill).  Keying the graph on quantity TYPES ([`VariableKey`])
# makes that structurally impossible — `ThermalEntropy` is not `Thermopower`, so the
# step never fires.  Operates over the type-keyed relations (those with identity
# slots); the supplied (untyped) slots are provided as `extras`, exactly as for the
# bag verbs.  See docs/design/type-keyed-interface.md §6.

"""
    TypedStep(relation, output::VariableKey, inputs::Vector{VariableKey})

One directed edge of the TYPE-keyed derivation graph: `relation` computes the
identity variable `output` (a [`VariableKey`](@ref)) from its other identity slots
`inputs` (plus its supplied slots, provided as `extras`), via the type-keyed
[`solve`](@ref).
"""
struct TypedStep
    relation::AbstractRelation
    output::VariableKey
    inputs::Vector{VariableKey}
end

const _TYPED_DERIV_STEPS = Ref{Union{Nothing,Vector{TypedStep}}}(nothing)

"""
    typed_derivation_steps() -> Vector{TypedStep}

Every candidate directed edge of the type-keyed derivation graph: for each
type-keyed relation and each of its identity slots, the edge computing that slot's
TYPE from the other identity slots.  Inequalities and symbol-only relations are
skipped (an inequality gives a saturation bound, not a derivation; a symbol-only
relation has no typed slots).  Built once and cached.
"""
function typed_derivation_steps()
    cached = _TYPED_DERIV_STEPS[]
    cached === nothing || return cached
    steps = TypedStep[]
    for rel in all_relations()
        rel isa AbstractInequality && continue
        idslots = [(s, k) for (s, k) in variable_slots(rel) if k !== nothing]
        isempty(idslots) && continue
        for (osym, okey) in idslots
            ins = [VariableKey(k) for (s, k) in idslots if s !== osym]
            push!(steps, TypedStep(rel, VariableKey(okey), ins))
        end
    end
    _TYPED_DERIV_STEPS[] = steps
    return steps
end

# Is a type available in `known` — aliasing InverseTemperature ⇄ Temperature (so a
# derived T is recognized as "known" when β was supplied, and vice versa, and neither
# is ever re-derived as a duplicate literal key that would trip `_slot_value`'s
# "both present" guard mid-chain).
_known(@nospecialize(ty::Type), known::Bag) = _slot_value(ty, known) !== nothing

# Reject a bag that supplies BOTH β and T up front (mirrors the bag verbs' guard),
# so the aliasing-aware checks below never encounter both literal keys at once.
function _check_one_temperature(bag::Bag)
    haskey(bag, VariableKey(InverseTemperature)) &&
        haskey(bag, VariableKey(Temperature)) &&
        error("bag has both InverseTemperature and Temperature — pass exactly one")
    return nothing
end

# Fire a typed step against the known bag + `extras`; the solved value or `nothing`
# (an input type is missing, or the relation is non-affine in the output / a
# supplied slot is absent — the real `solve` decides, so nothing is faked).  The
# whole attempt is guarded, so an unexpectedly-erroring step is skipped, never fatal.
function _try_typed_step(step::TypedStep, known::Bag, extras)
    try
        all(k -> _known(k.type, known), step.inputs) || return nothing
        return solve(step.relation, step.output.type, known; extras...)
    catch
        return nothing
    end
end

# Forward-chaining closure over VariableKey nodes: fire any step whose inputs are all
# known until nothing new is produced (stopping early once `stop` is known).  All
# "known" tests are aliasing-aware, so `known` never accumulates both β and T.
function _typed_chain!(known::Bag, extras, stop::Union{VariableKey,Nothing})
    used = TypedStep[]
    progress = true
    while progress && !(stop !== nothing && _known(stop.type, known))
        progress = false
        for step in typed_derivation_steps()
            _known(step.output.type, known) && continue
            v = _try_typed_step(step, known, extras)
            v === nothing && continue
            known[step.output] = v
            push!(used, step)
            progress = true
            stop !== nothing && _known(stop.type, known) && break
        end
    end
    return used
end

"""
    derivable(bag::Bag; extras...) -> Set{VariableKey}

Type-keyed [`derivable`](@ref): the quantity/field TYPES computable from the bag
`bag` (plus `extras` for supplied slots), through chains of the type-keyed
[`solve`](@ref).  Matched by TYPE, so no formula-symbol collision — a
`ThermalEntropy` in the bag can never pose as a `Thermopower`.
"""
function derivable(bag::Bag; extras...)
    known = copy(bag)
    _check_one_temperature(known)
    _typed_chain!(known, values(extras), nothing)
    return Set(keys(known))
end

"""
    derive(Q::Type, bag::Bag; extras...) -> value

Type-keyed [`derive`](@ref): the value of quantity/field TYPE `Q` computed from
`bag` (+ `extras`), or an error if unreachable.  Collision-proof — a relation fires
only when the ACTUAL quantity types it needs are present:

```julia
derive(PeltierCoefficient, bag(Thermopower => s, Temperature => t))     # t·s (Kelvin)
derive(PeltierCoefficient, bag(ThermalEntropy => s, Temperature => t))  # ERROR: unreachable
#   — KelvinRelation needs a Thermopower, not the entropy; the silent
#   entropy-as-Seebeck derivation the symbol-keyed graph allowed is impossible here.
```
"""
function derive(@nospecialize(Q::Type), bag::Bag; extras...)
    known = copy(bag)
    _check_one_temperature(known)
    v = _slot_value(Q, known)                         # aliasing-aware "already known"
    v === nothing || return something(v)
    _typed_chain!(known, values(extras), VariableKey(Q))
    v = _slot_value(Q, known)
    v === nothing && error(
        "derive: $(nameof(Q)) is not reachable from the bag by type-keyed relations " *
        "(a needed quantity type or supplied slot is absent).",
    )
    return something(v)
end

"""
    typed_derivation_graph() -> KnowledgeGraph{VariableKey}

The type-keyed derivation graph: one directed edge `input →[relation] output` per
(typed step, input), nodes are [`VariableKey`](@ref)s — the collision-proof
counterpart of [`derivation_graph`](@ref).  Structural (over-approximates, like its
symbol sibling; use [`derive`](@ref)`(Q, bag)` for honest reachability).
"""
function typed_derivation_graph()
    edges = TypedEdge{VariableKey}[]
    for step in typed_derivation_steps()
        rname = Symbol(nameof(typeof(step.relation)))
        for inp in step.inputs
            push!(edges, TypedEdge(rname, inp, step.output, string(rname), true))
        end
    end
    return KnowledgeGraph(edges)
end

# ─── Quantity-first navigation (a veneer over the typed graph) ───────────

"""
    reachable_quantities(q) -> Vector{Type}

The physical-quantity TYPES structurally reachable from `q` through the type-keyed
derivation graph — the quantity-first navigation over [`typed_derivation_graph`](@ref).
Accepts a quantity instance or its (concrete) type; the result is family-erased
(`Susceptibility{I}` → `Susceptibility`), de-duplicated, name-sorted, and includes
`q`'s own family (trivially reachable).

```julia
reachable_quantities(PartitionFunction())   # ⊇ [FreeEnergy, PartitionFunction] — F = −β⁻¹ ln Z
```

The dual of [`relations_constraining`](@ref) (a quantity → the laws it obeys); this
is a quantity → the other quantities its laws connect it to.

!!! warning "Structural, not computational"
    Reachability over [`typed_derivation_graph`](@ref) OVER-approximates what is
    actually computable (a hyperedge fans out to one edge per input, so a single
    known input already "reaches" the output; non-affine outputs are edges too).
    For the honest "can I compute it from these values" use [`derivable`](@ref)`(bag)`
    / [`derive`](@ref)`(Q, bag)`, which require every input and call the real `solve`.
"""
reachable_quantities(q::AbstractQuantity) = reachable_quantities(typeof(q))
function reachable_quantities(@nospecialize(Q::Type))
    fam = _family(Q)
    self = fam <: AbstractQuantity ? Type[fam] : Type[]
    g = typed_derivation_graph()
    start = VariableKey(Q)
    start in graph_nodes(g) || return self          # not a node ⇒ only its own family
    reached = graph_reachable(g, start)
    qs = unique!(Type[_family(k.type) for k in reached if k.type <: AbstractQuantity])
    return sort!(qs; by=T -> string(nameof(T)))
end

export DerivationStep,
    DerivationTrace,
    derivation_steps,
    derivation_graph,
    derivable,
    derive,
    reachable_quantities,
    TypedStep,
    typed_derivation_steps,
    typed_derivation_graph
