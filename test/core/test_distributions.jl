# Distribution vocabulary: type tree, constructors, ensemble weights,
# and the ThermalAverage marker.

using AbstractQAtlas
using AbstractQAtlas: solve

@testset "distribution type tree + constructors" begin
    @test MicroCanonical(1.5) isa AbstractDistribution
    @test MicroCanonical(1.5).ΔE == 0.0
    @test MicroCanonical(1.5; ΔE=0.2).ΔE == 0.2
    @test Canonical(2.0) isa AbstractDistribution
    @test Canonical(2.0).β == 2.0
    @test Canonical(; T=0.5).β == 2.0            # β-or-T convention
    @test Canonical(; β=1//2).β === 1//2         # exact types preserved
    gc = GrandCanonical(; T=1.0, μ=0.3)
    @test gc isa AbstractDistribution
    @test gc.β == 1.0 && gc.μ == 0.3
    sq = Squeezed(0.8)
    @test sq isa AbstractDistribution
    @test sq.r == 0.8 && sq.φ == 0.0
    @test Squeezed(0.8; φ=π / 4).φ ≈ π / 4
    @test Fermionic() isa ParticleStatistics
    @test Bosonic() isa ParticleStatistics
end

@testset "ensemble weights" begin
    # canonical weights over a spectrum sum to Z — cross-checked through
    # the FreeEnergyFromZ relation (two independent code paths).
    levels = [0.0, 1.3]
    β = 0.9
    Z = sum(ensemble_weight(Canonical(β), E) for E in levels)
    @test Z ≈ 1 + exp(-β * 1.3) atol = 1e-14
    F = solve(FreeEnergyFromZ(), Val(:f); Z=Z, β=β)
    @test solve(FreeEnergyFromZ(), Val(:Z); f=F, β=β) ≈ Z rtol = 1e-12

    # grand canonical reduces to canonical at μ = 0, and shifts by e^{βμN}
    gc = GrandCanonical(β, 0.0)
    @test ensemble_weight(gc, 1.3; N=2) == ensemble_weight(Canonical(β), 1.3)
    gc2 = GrandCanonical(β, 0.5)
    @test ensemble_weight(gc2, 1.3; N=2) ≈ exp(β * 0.5 * 2) * exp(-β * 1.3) atol = 1e-14

    # microcanonical indicator: counts exactly the states in the window
    mc = MicroCanonical(1.0; ΔE=0.4)
    spectrum = [0.0, 0.85, 1.0, 1.19, 1.5]
    @test sum(ensemble_weight(mc, E) for E in spectrum) == 3   # 0.85, 1.0, 1.19
    @test ensemble_weight(mc, 1.21) == 0
end

@testset "ThermalAverage marker" begin
    ta = ThermalAverage(Magnetization(:z), Canonical(1.0))
    @test ta isa AbstractQuantity                 # composes with fetch
    @test ta.quantity === Magnetization(:z)
    @test ta.distribution.β == 1.0
    # tensor traits pass through the marker
    @test indices(ta) == (:z,)
    @test indices(ThermalAverage(Susceptibility(:x, :y), GrandCanonical(1.0, 0.0))) ==
        (:x, :y)
    @test tensor_rank(ThermalAverage(Susceptibility(:x, :y), Canonical(1.0))) == 2
    @test indices(ThermalAverage(SpecificHeat(), Canonical(1.0))) == ()
end
