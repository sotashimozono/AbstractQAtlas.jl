# Quantum statistics vs structural identities and limits — every check
# an independent expectation (exact symmetry, exact identity between the
# two statistics, classical limit, closed-form trig identities).

using AbstractQAtlas

@testset "Fermi–Dirac structure" begin
    β = 1.7
    μ = 0.4
    # exact half filling at the chemical potential
    @test occupation(Fermionic(), μ; β=β, μ=μ) == 0.5
    # particle–hole symmetry: n(μ+δ) + n(μ−δ) = 1 (exact identity)
    for δ in (0.1, 0.7, 3.0)
        s =
            occupation(Fermionic(), μ + δ; β=β, μ=μ) +
            occupation(Fermionic(), μ - δ; β=β, μ=μ)
        @test s ≈ 1.0 atol = 1e-14
    end
    # bounded in [0, 1]
    @test 0 < occupation(Fermionic(), -10.0; β=β, μ=μ) <= 1
    # T-form kwarg equals β-form
    @test occupation(Fermionic(), 1.0; T=1 / β, μ=μ) ==
        occupation(Fermionic(), 1.0; β=β, μ=μ)
end

@testset "Bose–Einstein structure" begin
    β = 1.7
    # domain guard: ε ≤ μ is rejected
    @test_throws ErrorException occupation(Bosonic(), 0.0; β=β, μ=0.0)
    @test_throws ErrorException occupation(Bosonic(), -0.5; β=β, μ=0.0)
    # divergence as ε → μ⁺ (monotone growth toward the condensation edge)
    @test occupation(Bosonic(), 0.01; β=β) > occupation(Bosonic(), 0.1; β=β)
    @test occupation(Bosonic(), 1e-4; β=β) > 1e3 * occupation(Bosonic(), 1.0; β=β)
end

@testset "exact structural identity: n_B(ε) − n_F(ε) = 2 n_B(2ε) at μ = 0" begin
    β = 0.9
    for ε in (0.2, 1.0, 3.7)
        lhs = occupation(Bosonic(), ε; β=β) - occupation(Fermionic(), ε; β=β)
        rhs = 2 * occupation(Bosonic(), 2ε; β=β)
        @test lhs ≈ rhs atol = 1e-14
    end
end

@testset "classical (Boltzmann) limit for β(ε−μ) ≫ 1" begin
    β = 1.0
    ε = 12.0   # β(ε−μ) = 12 ⇒ quantum corrections ~ e^{−12} ≈ 6e-6 relative
    nB = occupation(Bosonic(), ε; β=β)
    nF = occupation(Fermionic(), ε; β=β)
    ncl = boltzmann_occupation(ε; β=β)
    @test isapprox(nB, ncl; rtol=1e-5)
    @test isapprox(nF, ncl; rtol=1e-5)
    # and the quantum corrections bracket the classical value:
    @test nF < ncl < nB
end

@testset "squeezed-state moments" begin
    # r = 0 is the vacuum: Var = 1/2 each, zero photons
    v0 = squeezed_variances(0.0)
    @test v0.x == 0.5 && v0.p == 0.5
    @test squeezed_mean_photons(0.0) == 0.0
    for r in (0.3, 0.8, 2.0)
        v = squeezed_variances(r)
        # minimum-uncertainty saturation for every r (exact closed form)
        @test v.x * v.p ≈ 1 / 4 atol = 1e-14
        # squeezing direction: x squeezed below vacuum, p anti-squeezed
        @test v.x < 0.5 < v.p
        # independent route to ⟨n⟩: (Var(x)+Var(p))/2 − 1/2 == sinh²r
        # (⟨x²⟩+⟨p²⟩ = 2⟨n⟩+1 — a different identity than the definition)
        @test squeezed_mean_photons(r) ≈ (v.x + v.p) / 2 - 1 / 2 atol = 1e-12
    end
end
