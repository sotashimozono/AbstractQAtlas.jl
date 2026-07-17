# test/relations/test_region_entropy.jl — §5 Region set layer + region-keyed
# entanglement-entropy auto-discovery (design doc §5/§8b, Phase-2 P1).

using AbstractQAtlas
using Test
const AQ = AbstractQAtlas

@testset "Region set algebra (dimension-agnostic)" begin
    A, B = Region(1, 2), Region(3, 4)
    @test disjoint(A, B)
    @test !disjoint(A, Region(2, 3))
    @test A ∪ B == Region(1, 2, 3, 4)
    @test A ∩ Region(2, 3) == Region(2)
    @test Region(1) ⊆ A && !(Region(1, 5) ⊆ A)
    @test length(A) == 2 && !isempty(A) && isempty(Region{Int}(Set{Int}()))
    # value-based equality/hash — order-insensitive
    @test Region(2, 1) == Region(1, 2)
    @test hash(Region(2, 1)) == hash(Region(1, 2))
    # ND sites: 2D tuples work identically (the set layer is ND from day one)
    a2, b2 = Region((1, 1), (1, 2)), Region((2, 1))
    @test disjoint(a2, b2)
    @test a2 ∪ b2 == Region((1, 1), (1, 2), (2, 1))
end

@testset "entanglement_entropy() region-keyed bag" begin
    @test entanglement_entropy(1, 2) == entanglement_entropy(Region(1, 2))                 # value-based key
    @test entanglement_entropy(1, 2) isa VariableKey
    @test entanglement_entropy(1, 2).type === VonNeumannEntropy
    @test entanglement_entropy(1, 2).support isa RegionSupport
    @test entanglement_entropy(2, 1) == entanglement_entropy(1, 2)                          # region order-insensitive
    @test bag(entanglement_entropy(1) => 0.5, entanglement_entropy(1, 2) => 1.0)[entanglement_entropy(
        1
    )] == 0.5
end

@testset "region_report: auto-discovery of subadditivity + Araki–Lieb" begin
    # subadditive bag: S(A)+S(B) ≥ S(A∪B) and S(A∪B) ≥ |S(A)−S(B)| — all hold
    good = bag(
        entanglement_entropy(1) => 0.7,
        entanglement_entropy(2) => 0.7,
        entanglement_entropy(1, 2) => 1.0,
    )
    rep = region_report(good)
    @test length(rep) == 2                                       # Subadditivity + ArakiLieb, one region pair
    @test all(r -> r.pass, rep)
    @test region_check_all(good)
    @test Set(nameof(typeof(r.relation)) for r in rep) == Set((:Subadditivity, :ArakiLieb))

    # a NEGATIVE mutual information (unphysical — a broken calc) is caught, on the
    # right region pair, with no A/B/AB hand-labeling
    bad = bag(
        entanglement_entropy(1) => 0.5,
        entanglement_entropy(2) => 0.5,
        entanglement_entropy(1, 2) => 1.5,
    )   # I(A:B) = −0.5
    @test !region_check_all(bad)
    viol = [r for r in region_report(bad) if !r.pass]
    @test length(viol) == 1
    @test viol[1].relation isa Subadditivity
    @test Set(viol[1].regions) == Set((Region(1), Region(2)))
    @test viol[1].slack ≈ -0.5

    # 3 regions ⇒ every disjoint pair whose union is present is auto-discovered
    b3 = bag(
        entanglement_entropy(1) => 0.6,
        entanglement_entropy(2) => 0.6,
        entanglement_entropy(3) => 0.6,
        entanglement_entropy(1, 2) => 1.0,
        entanglement_entropy(1, 3) => 1.0,
        entanglement_entropy(2, 3) => 1.0,
    )
    @test length(region_report(b3)) == 6                         # 3 disjoint pairs × 2 relations
    @test region_check_all(b3)

    # empty match ⇒ false, never a silent green
    @test !region_check_all(bag(entanglement_entropy(1) => 0.5))
    @test isempty(region_report(bag(entanglement_entropy(1) => 0.5)))
    # a non-disjoint pair (or a missing union / missing S(B)) yields no instance
    @test isempty(
        region_report(
            bag(entanglement_entropy(1) => 0.5, entanglement_entropy(1, 2) => 1.0)
        ),
    )
end

@testset "region_report: strong subadditivity (triples) + mutual_information" begin
    ee = entanglement_entropy
    # SSA-valid bag (atoms A={1}, B={2}, C={3}): S(A∪B)+S(B∪C) ≥ S(A∪B∪C)+S(B)
    b = bag(
        ee(1) => 0.5,
        ee(2) => 0.5,
        ee(3) => 0.5,
        ee(1, 2) => 1.0,
        ee(2, 3) => 1.0,
        ee(1, 2, 3) => 1.2,
    )
    ssa = [r for r in region_report(b) if r.relation isa StrongSubadditivity]
    @test length(ssa) == 1                                # one triple, B={2} the shared middle
    @test ssa[1].pass && ssa[1].slack ≈ 0.3
    @test Set(ssa[1].regions) == Set((Region(1), Region(2), Region(3)))
    @test region_check_all(b)                             # SSA + subadditivity + Araki–Lieb all hold

    # a FULLY-connected triple: all 3 choices of the shared middle B are valid, so the
    # enumeration must emit exactly 3 distinct SSA instances (no dup, no omission)
    full = bag(
        ee(1) => 0.5,
        ee(2) => 0.5,
        ee(3) => 0.5,
        ee(1, 2) => 1.0,
        ee(1, 3) => 1.0,
        ee(2, 3) => 1.0,
        ee(1, 2, 3) => 1.2,
    )
    ssa_full = [r for r in region_report(full) if r.relation isa StrongSubadditivity]
    @test length(ssa_full) == 3                           # one per choice of the middle B
    @test length(unique(r.regions[2] for r in ssa_full)) == 3      # 3 distinct middles
    @test all(r -> r.pass, ssa_full)

    # a broken calc — negative conditional mutual information I(A:C|B) < 0 — is caught
    bad = bag(
        ee(1) => 0.5,
        ee(2) => 0.5,
        ee(3) => 0.5,
        ee(1, 2) => 1.0,
        ee(2, 3) => 1.0,
        ee(1, 2, 3) => 1.8,
    )
    @test !region_check_all(bad)
    viol = [r for r in region_report(bad) if r.relation isa StrongSubadditivity && !r.pass]
    @test length(viol) == 1 && viol[1].slack ≈ -0.3

    # mutual_information helper: I(A:B) = S(A) + S(B) − S(A∪B)
    @test mutual_information(b, Region(1), Region(2)) ≈ 0.0                 # 0.5 + 0.5 − 1.0
    @test mutual_information(
        bag(ee(1) => 0.7, ee(2) => 0.7, ee(1, 2) => 1.0), Region(1), Region(2)
    ) ≈ 0.4
    @test_throws ErrorException mutual_information(bag(ee(1) => 0.5), Region(1), Region(2))
end

@testset "multipartite region helpers: CMI, tripartite info, KP topological EE" begin
    ee = entanglement_entropy
    A, B, C = Region(1), Region(2), Region(3)
    # pairwise-DISTINCT singles and pairs so a positional swap among {S_A,S_B,S_C} or
    # {S_AB,S_AC,S_BC} — the failure mode of the generator destructuring — cannot hide
    b = bag(
        ee(1) => 0.5,
        ee(2) => 0.6,
        ee(3) => 0.7,
        ee(1, 2) => 1.0,
        ee(1, 3) => 1.1,
        ee(2, 3) => 1.2,
        ee(1, 2, 3) => 1.2,
    )
    # I(A:C|B) = S(A∪B)+S(B∪C)−S(A∪B∪C)−S(B) = 1.0+1.2−1.2−0.6   (the SSA slack; S(B)=S₂
    # is distinct from S₁/S₃, so picking the wrong middle region would change the answer)
    cmi = conditional_mutual_information(b, A, B, C)
    @test cmi ≈ 0.4
    # I₃ = S_A+S_B+S_C−S_AB−S_AC−S_BC+S_ABC = 1.8−3.3+1.2 ; Kitaev–Preskill γ = −I₃
    @test tripartite_information(b, A, B, C) ≈ -0.3
    @test topological_entanglement_entropy(b, A, B, C) ≈ 0.3
    @test topological_entanglement_entropy(b, A, B, C) ≈ -tripartite_information(b, A, B, C)
    # CMI agrees with the auto-discovered StrongSubadditivity slack for the same triple
    @test any(r -> r.relation isa StrongSubadditivity && r.slack ≈ cmi, region_report(b))
    # a missing entropy is a loud error, never silently wrong
    @test_throws ErrorException tripartite_information(bag(ee(1) => 0.5), A, B, C)
    @test_throws ErrorException conditional_mutual_information(
        bag(ee(1, 2) => 1.0), A, B, C
    )
    @test_throws ErrorException topological_entanglement_entropy(bag(ee(1) => 0.5), A, B, C)
    # overlapping regions are NOT a tripartition (unions collapse ⇒ a physical-looking
    # but meaningless number) — caught by the disjointness guard, on a bag where every
    # needed entropy IS present so only the guard can be firing
    @test_throws ErrorException conditional_mutual_information(b, A, B, A)   # C == A
    @test_throws ErrorException tripartite_information(b, A, B, A)
    @test_throws ErrorException topological_entanglement_entropy(b, A, B, A)
end
