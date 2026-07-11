# Scaling relations vs INDEPENDENT exact exponent sets.
#
# The 2D Ising rationals and the mean-field set are exact, independently
# known values (Onsager/Yang lattice solutions; Landau theory) — the
# relations must hold with residual ≡ 0 in EXACT arithmetic, which also
# pins the no-float-promotion contract.

using AbstractQAtlas
using AbstractQAtlas: residual, check, solve

const ISING2D = (α=0//1, β=1//8, γ=7//4, δ=15//1, ν=1//1, η=1//4)
const MEANFIELD = (α=0//1, β=1//2, γ=1//1, δ=3//1, ν=1//2, η=0//1)

@testset "2D Ising exact rationals: residuals ≡ 0//1, types preserved" begin
    r1 = residual(Rushbrooke(); α=ISING2D.α, β=ISING2D.β, γ=ISING2D.γ)
    r2 = residual(Widom(); β=ISING2D.β, γ=ISING2D.γ, δ=ISING2D.δ)
    r3 = residual(Fisher(); γ=ISING2D.γ, ν=ISING2D.ν, η=ISING2D.η)
    r4 = residual(Josephson(); α=ISING2D.α, ν=ISING2D.ν, d=2)
    for r in (r1, r2, r3, r4)
        @test r isa Rational       # exact-arithmetic contract
        @test r == 0//1            # exact, not ≈
    end
    @test exponents_consistent(ISING2D; d=2)          # atol = 0: exact gate
    @test all(iszero, values(exponent_residuals(ISING2D; d=2)))
end

@testset "mean-field set at the upper critical dimension d = 4" begin
    @test exponents_consistent(MEANFIELD; d=4)
    # hyperscaling FAILS off the upper critical dimension — that failure
    # is physics (mean-field violates Josephson for d < 4), so the gate
    # must flag it:
    @test !check(Josephson(); α=MEANFIELD.α, ν=MEANFIELD.ν, d=3)
end

@testset "3D Ising bootstrap values within quoted precision" begin
    # Kos–Poland–Simmons-Duffin–Vichi (2016) determine (Δ_σ, Δ_ε); the
    # standard exponent set derives from them, so this is a consistency
    # check of our arithmetic against the published rounded values (the
    # genuinely independent exactness test is the 2D rational case above).
    nt3 = (α=0.11009, β=0.32642, γ=1.23708, δ=4.78984, ν=0.62999, η=0.03631)
    # residuals limited by the 5-digit rounding of the published values:
    @test exponents_consistent(nt3; d=3, atol=2e-4)
end

@testset "solve: every solvable variable round-trips exactly" begin
    @test solve(Rushbrooke(), Val(:α); β=1//8, γ=7//4) == 0//1
    @test solve(Rushbrooke(), Val(:β); α=0//1, γ=7//4) == 1//8
    @test solve(Rushbrooke(), Val(:γ); α=0//1, β=1//8) == 7//4
    @test solve(Widom(), Val(:β); γ=7//4, δ=15//1) == 1//8
    @test solve(Widom(), Val(:γ); β=1//8, δ=15//1) == 7//4
    @test solve(Widom(), Val(:δ); β=1//8, γ=7//4) == 15//1
    @test solve(Fisher(), Val(:γ); ν=1//1, η=1//4) == 7//4
    @test solve(Fisher(), Val(:ν); γ=7//4, η=1//4) == 1//1
    @test solve(Fisher(), Val(:η); γ=7//4, ν=1//1) == 1//4
    @test solve(Josephson(), Val(:α); ν=1//1, d=2) == 0//1
    @test solve(Josephson(), Val(:ν); α=0//1, d=2) == 1//1
    @test solve(Josephson(), Val(:d); α=0//1, ν=1//1) == 2//1
    # solve-then-residual ≡ 0 (definitional round trip)
    γ = solve(Widom(), Val(:γ); β=1//8, δ=15//1)
    @test residual(Widom(); β=1//8, γ=γ, δ=15//1) == 0//1
end

@testset "a wrong exponent set fails the gate" begin
    bad = (α=0//1, β=1//8, γ=3//2, δ=15//1, ν=1//1, η=1//4)   # γ ≠ 7/4
    @test !exponents_consistent(bad; d=2)
    @test !check(Widom(); β=bad.β, γ=bad.γ, δ=bad.δ)
end
