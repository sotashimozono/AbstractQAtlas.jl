# relations/region_entropy.jl — auto-discovery of the entanglement-entropy
# inequalities over the REGIONS present in a bag (design §5/§8b, Phase-2 P1).
#
# The bipartite entropy inequalities (subadditivity, Araki–Lieb) hold for ANY
# disjoint pair of regions (A, B).  Keyed on a Region support (`entropy(A)`), they
# become auto-discoverable: `region_report` scans a bag of region-entropies, finds
# every disjoint (A, B) whose S(A), S(B), S(A∪B) are all present, and checks the
# inequalities on that instance — no hand-labeled A/B/AB.  The relations' scalar
# kernels (Subadditivity, ArakiLieb) are reused verbatim; this is the region-matching
# layer over them.  SSA / multipartite / §8b index-unification follow in P2 on this
# same matcher.

"""
    RegionReportRow

One row of a [`region_report`](@ref): the `relation` (an entropy inequality), the
disjoint region pair `(A, B)` it was auto-instantiated on, the `slack` (its
[`residual`](@ref); `≥ 0` ⇔ satisfied), and `pass`.
"""
struct RegionReportRow
    relation::AbstractRelation
    regions::Tuple{Region,Region}
    slack::Number
    pass::Bool
end
export RegionReportRow

# the von Neumann region-entropies present in a bag: Region → S(Region)
function _region_entropies(b::Bag)
    return Dict(
        k.support.region => v for
        (k, v) in b if k.type === VonNeumannEntropy && k.support isa RegionSupport
    )
end

"""
    region_report(b::Bag; atol=0) -> Vector{RegionReportRow}

Auto-discover the bipartite entanglement-entropy inequalities over the REGIONS in a
bag of region-keyed entropies (`bag(entropy(A) => s_A, …)`): for every disjoint pair
(A, B) whose `S(A)`, `S(B)`, `S(A∪B)` are all present, check **subadditivity**
`I(A:B) = S(A) + S(B) − S(A∪B) ≥ 0` and the **Araki–Lieb** bound
`S(A∪B) ≥ |S(A) − S(B)|` — one row per (relation, unordered region pair).  The
region twin of [`relation_report`](@ref): a negative mutual information (a broken
MPS/ED entanglement calculation) is caught for whichever regions expose it, with no
A/B/AB hand-labeling.

```julia
b = bag(entropy(1) => 0.7, entropy(2) => 0.7, entropy(1, 2) => 1.0)   # S(A), S(B), S(A∪B)
all(row -> row.pass, region_report(b))     # true — S is subadditive here
```
"""
function region_report(b::Bag; atol=0)
    ents = _region_entropies(b)
    regions = collect(keys(ents))
    out = RegionReportRow[]
    for i in eachindex(regions), j in (i + 1):lastindex(regions)
        A, B = regions[i], regions[j]
        disjoint(A, B) || continue
        haskey(ents, A ∪ B) || continue
        S_A, S_B, S_AB = ents[A], ents[B], ents[A ∪ B]
        for rel in (Subadditivity(), ArakiLieb())
            s = residual(rel; S_A=S_A, S_B=S_B, S_AB=S_AB)
            push!(out, RegionReportRow(rel, (A, B), s, _passes(rel, s, atol)))
        end
    end
    return out
end
export region_report

"""
    region_check_all(b::Bag; atol=0) -> Bool

`true` iff every bipartite entropy inequality auto-discovered by
[`region_report`](@ref) holds on the bag `b` — and at least one instance was found
(an empty match is `false`, never a silent green).
"""
function region_check_all(b::Bag; atol=0)
    rep = region_report(b; atol=atol)
    return !isempty(rep) && all(row -> row.pass, rep)
end
export region_check_all
