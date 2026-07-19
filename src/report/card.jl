# report/card.jl — the report / verification-CARD contract (the ecosystem seam).
#
# Where `fetch(model, quantity, bc)` RETRIEVES a value, `report(...)` PACKAGES a
# computed one — an oracle's measurement plus its provenance and error bar — into a
# schema-valid [`Card`](@ref).  A reporter package (ClassicalMonteCarlo, WilsonNRG,
# MPS-QMC, …) depends on AbstractQAtlas ONLY and emits cards; the registry is a data
# sink, an analyzer reconciles + gates, a documenter renders.  Reporters push DATA,
# not a code dependency — so the ecosystem graph stays cycle-free.
#
# The card SCHEMA (v2) is shared with QAtlas's black-box `verify(...)` cards
# (QAtlas test/util/verify.jl): hub / route / mechanism / independence / status /
# subject / independent[] / atol / refs, hardened here with an explicit `provenance`
# discriminant and an `error_bar`, and NaN/Inf-safe (a non-finite subject becomes
# `status = :divergent`, never a raw NaN token).  One schema ⇒ one consumer renders
# verify-cards and report-cards uniformly.

const CARD_SCHEMA_VERSION = 2

"""
    REPORT_ROUTES

The recognized `route`s a [`Card`](@ref) may carry — the *how* of the value.  The
cross-check routes are adopted verbatim from QAtlas's schema-v2 `verify`
vocabulary; the measurement routes (`:monte_carlo`, `:dmrg`, `:mps_qmc`, `:nrg`,
`:tpq`) are for the oracle reporters that push cards into the registry.
"""
const REPORT_ROUTES = (
    # cross-check routes — a value corroborated against an independent one
    :ed_finite_size,
    :sum_rule,
    :delegation_invariant,
    :limiting_case,
    :literature_value,
    :second_closed_form,
    # measurement routes — an oracle reporting its own computed value
    :monte_carlo,
    :dmrg,
    :mps_qmc,
    :nrg,
    :tpq,
)

# Only these routes are mechanically INDEPENDENT of an atlas's closed form, so a
# corroboration by one is "structural"; any other cross-check is "asserted" (QAtlas
# review B1).  A bare measurement (no `independent`) is neither — it is `:measured`.
const _CARD_STRUCTURAL_ROUTES = (:ed_finite_size, :second_closed_form, :literature_value)

"""
    Card

One schema-v2 verification/report card: a computed `subject` value for
`hub = "TypeName(model)/TypeName(quantity)/TypeName(bc)"`, obtained via `route`,
with its `error_bar`, `independence` class, `status`, any `independent` cross-check
values, `atol`, `refs`, and a `provenance` discriminant.  Build one with
[`report`](@ref); serialize a stream with [`card_jsonl`](@ref).

`status` is `:divergent` (and `subject` is `nothing`) when the reported value is
non-finite — a NaN/Inf is never emitted as a raw token.
"""
struct Card
    schema_version::Int
    hub::String
    route::Symbol
    mechanism::String
    independence::Symbol                 # :structural | :asserted | :measured
    status::Symbol                       # :ok | :divergent
    subject::Union{Float64,Nothing}
    error_bar::Union{Float64,Nothing}
    independent::Vector{Float64}
    atol::Float64
    refs::Vector{String}
    provenance::String
end

function Base.show(io::IO, c::Card)
    val = c.subject === nothing ? "divergent" : string(c.subject)
    eb = c.error_bar === nothing ? "" : " ± " * string(c.error_bar)
    print(io, "Card(", c.hub, " = ", val, eb, "  [", c.route, "/", c.independence, "])")
    return nothing
end

# The type NAME of a model/quantity/bc, whether given as an instance or a type.
_card_typename(@nospecialize(x)) = string(nameof(x isa Type ? x : typeof(x)))

function _card_hub(model, quantity, bc)
    return string(
        _card_typename(model), "/", _card_typename(quantity), "/", _card_typename(bc)
    )
end

# A real `Float64` if finite; `nothing` if NaN/Inf; a Complex with a negligible
# imaginary part is taken as its real part (correlators of Hermitian operators carry
# a few ULPs of round-off) — mirrors QAtlas's JSON-safe numeric handling.
function _card_finite(x)
    r = if x isa Complex
        abs(imag(x)) <= 1e-9 * max(1.0, abs(real(x))) ? real(x) : NaN
    else
        x
    end
    v = float(r)
    return isfinite(v) ? Float64(v) : nothing
end

# The independence class from the route + whether a cross-check was supplied: a bare
# measurement is `:measured`; a corroborated value is `:structural` iff its route is
# mechanically independent, else `:asserted`.
function _card_independence(route::Symbol, independent)
    isempty(independent) && return :measured
    return route in _CARD_STRUCTURAL_ROUTES ? :structural : :asserted
end

"""
    report(model, quantity, bc;
           value, route, provenance,
           err=nothing, mechanism="", independent=(), atol=0, refs=()) -> Card

Package a computed `value` for the `(model, quantity, bc)` triple into a
schema-valid [`Card`](@ref) — the reporter-facing sibling of
[`fetch`](@ref AbstractQAtlas.fetch) (fetch RETRIEVES a value; report PACKAGES one).
The card's `hub` is `"TypeName(model)/TypeName(quantity)/TypeName(bc)"` (instances or
types are both accepted).  `route` must be one of [`REPORT_ROUTES`](@ref);
`provenance` names what produced the value (the reporter package + method).  A
non-finite `value` yields `status = :divergent` with `subject = nothing`, never a
raw NaN.

`value` (and each `independent` entry) is a scalar `Real`/`Complex`; a complex value
with a negligible imaginary part is taken as real.

```julia
report(TFIM(1.0, 0.5), VonNeumannEntropy(), PBC(64);
       value = 0.87, err = 0.01, route = :monte_carlo,
       provenance = "ClassicalMonteCarlo.jl@metropolis", refs = ["Calabrese2004"])
```
"""
function report(
    model,
    quantity,
    bc;
    value,
    route::Symbol,
    provenance::AbstractString,
    err::Union{Nothing,Real}=nothing,
    mechanism::AbstractString="",
    independent=(),
    atol::Real=0,
    refs=(),
)
    route in REPORT_ROUTES ||
        error("report: route must be one of $(REPORT_ROUTES); got $(repr(route))")
    subject = _card_finite(value)
    status = subject === nothing ? :divergent : :ok
    ind = Float64[]
    for x in independent
        f = _card_finite(x)
        f === nothing || push!(ind, f)     # a non-finite cross-check corroborates nothing
    end
    return Card(
        CARD_SCHEMA_VERSION,
        _card_hub(model, quantity, bc),
        route,
        String(mechanism),
        _card_independence(route, ind),          # classify by the SURVIVING cross-checks
        status,
        subject,
        err === nothing ? nothing : _card_finite(err),
        ind,
        Float64(atol),
        String[string(r) for r in refs],
        String(provenance),
    )
end

# NaN/Inf-safe JSON number: a finite value as itself, `nothing`/non-finite as null.
_card_json_num(::Nothing) = "null"
_card_json_num(x::Real) = isfinite(x) ? string(float(x)) : "null"
_card_json_arr(xs) = "[" * join((_card_json_num(x) for x in xs), ",") * "]"
_card_json_strarr(xs) = "[" * join((_json_str(x) for x in xs), ",") * "]"

"""
    card_jsonl([io=stdout], cards) -> nothing

Stream `cards` (any iterable of [`Card`](@ref)s) as JSONL — one schema-v2 object per
line — for the registry sink / a documenter.  NaN/Inf-safe (a non-finite
`subject`/`error_bar` is `null`, the card's `status` already `:divergent`).  Matches
the ecosystem `*_jsonl` idiom ([`graph_jsonl`](@ref)) so one consumer renders
models ⊕ quantities ⊕ derivations ⊕ cards.
"""
function card_jsonl(io::IO, cards)
    for c in cards
        print(
            io,
            "{\"schema_version\":", c.schema_version,
            ",\"hub\":", _json_str(c.hub),
            ",\"route\":", _json_str(c.route),
            ",\"mechanism\":", _json_str(c.mechanism),
            ",\"independence\":", _json_str(c.independence),
            ",\"status\":", _json_str(c.status),
            ",\"subject\":", _card_json_num(c.subject),
            ",\"error_bar\":", _card_json_num(c.error_bar),
            ",\"independent\":", _card_json_arr(c.independent),
            ",\"atol\":", _card_json_num(c.atol),
            ",\"refs\":", _card_json_strarr(c.refs),
            ",\"provenance\":", _json_str(c.provenance),
            "}\n",
        )
    end
    return nothing
end
card_jsonl(cards) = card_jsonl(stdout, cards)

export Card, report, card_jsonl, REPORT_ROUTES
