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

@testset "entropy() region-keyed bag" begin
    @test entropy(1, 2) == entropy(Region(1, 2))                 # value-based key
    @test entropy(1, 2) isa VariableKey
    @test entropy(1, 2).type === VonNeumannEntropy
    @test entropy(1, 2).support isa RegionSupport
    @test entropy(2, 1) == entropy(1, 2)                          # region order-insensitive
    @test bag(entropy(1) => 0.5, entropy(1, 2) => 1.0)[entropy(1)] == 0.5
end

@testset "region_report: auto-discovery of subadditivity + Araki–Lieb" begin
    # subadditive bag: S(A)+S(B) ≥ S(A∪B) and S(A∪B) ≥ |S(A)−S(B)| — all hold
    good = bag(entropy(1) => 0.7, entropy(2) => 0.7, entropy(1, 2) => 1.0)
    rep = region_report(good)
    @test length(rep) == 2                                       # Subadditivity + ArakiLieb, one region pair
    @test all(r -> r.pass, rep)
    @test region_check_all(good)
    @test Set(nameof(typeof(r.relation)) for r in rep) == Set((:Subadditivity, :ArakiLieb))

    # a NEGATIVE mutual information (unphysical — a broken calc) is caught, on the
    # right region pair, with no A/B/AB hand-labeling
    bad = bag(entropy(1) => 0.5, entropy(2) => 0.5, entropy(1, 2) => 1.5)   # I(A:B) = −0.5
    @test !region_check_all(bad)
    viol = [r for r in region_report(bad) if !r.pass]
    @test length(viol) == 1
    @test viol[1].relation isa Subadditivity
    @test Set(viol[1].regions) == Set((Region(1), Region(2)))
    @test viol[1].slack ≈ -0.5

    # 3 regions ⇒ every disjoint pair whose union is present is auto-discovered
    b3 = bag(
        entropy(1) => 0.6,
        entropy(2) => 0.6,
        entropy(3) => 0.6,
        entropy(1, 2) => 1.0,
        entropy(1, 3) => 1.0,
        entropy(2, 3) => 1.0,
    )
    @test length(region_report(b3)) == 6                         # 3 disjoint pairs × 2 relations
    @test region_check_all(b3)

    # empty match ⇒ false, never a silent green
    @test !region_check_all(bag(entropy(1) => 0.5))
    @test isempty(region_report(bag(entropy(1) => 0.5)))
    # a non-disjoint pair (or a missing union / missing S(B)) yields no instance
    @test isempty(region_report(bag(entropy(1) => 0.5, entropy(1, 2) => 1.0)))
end
