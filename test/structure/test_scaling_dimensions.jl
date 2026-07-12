# The RG-eigenvalue origin of the critical exponents: `critical_exponents`
# derives (α,β,γ,δ,ν,η) from (y_t, y_h, d), and the four scaling laws in
# relations/scaling.jl become CONSEQUENCES — their residual is identically 0
# for the derived set, at any eigenvalues.  The test checks the derivation
# against independent expectations (the scaling @relations themselves, known
# exponent tables, the inverse map), NOT against its own output.

using AbstractQAtlas
using AbstractQAtlas:
    ScalingDimensions,
    critical_exponents,
    critical_exponent,
    scaling_dimensions,
    residual,
    Rushbrooke,
    Widom,
    Fisher,
    Josephson,
    exponents_consistent

# a spread of EXACT rational RG data (various y_t, y_h, d), including the two
# physical fixed points 2D-Ising (1, 15//8, 2) and 3D-percolation-like sets
const _RG_SETS = (
    (1 // 1, 15 // 8, 2),      # 2D Ising (exact)
    (3 // 2, 7 // 4, 3),
    (2 // 1, 9 // 5, 3),
    (5 // 4, 11 // 6, 2),
    (4 // 3, 21 // 11, 3),
)

@testset "the four scaling laws are IDENTITIES in (y_t, y_h, d)" begin
    # feed the DERIVED exponents to the independent scaling @relations — every
    # residual must be exactly 0//1, for every eigenvalue set.  This is the
    # structural claim: Rushbrooke/Widom/Fisher/Josephson are not axioms but
    # consequences of two-eigenvalue homogeneity + hyperscaling.
    for (yt, yh, dd) in _RG_SETS
        e = critical_exponents(ScalingDimensions(yt, yh, dd))
        @test residual(Rushbrooke(); α=e.α, β=e.β, γ=e.γ) == 0 // 1
        @test residual(Widom(); β=e.β, γ=e.γ, δ=e.δ) == 0 // 1
        @test residual(Fisher(); γ=e.γ, ν=e.ν, η=e.η) == 0 // 1
        @test residual(Josephson(); α=e.α, ν=e.ν, d=dd) == 0 // 1
        # and the domain-wide gate agrees
        @test exponents_consistent(e; d=dd)
        # exactness: Rational eigenvalues ⇒ Rational exponents
        @test all(v -> v isa Rational, values(e))
    end
end

@testset "2D Ising exponents fall out of (1, 15//8, 2)" begin
    e = critical_exponents(ScalingDimensions(1 // 1, 15 // 8, 2))
    @test e == (α=0 // 1, β=1 // 8, γ=7 // 4, δ=15 // 1, ν=1 // 1, η=1 // 4)
    # single-exponent accessor agrees with the full set
    for name in (:α, :β, :γ, :δ, :ν, :η)
        @test critical_exponent(name, ScalingDimensions(1 // 1, 15 // 8, 2)) == e[name]
    end
end

@testset "mean-field / Gaussian eigenvalues at the upper critical dimension" begin
    # Gaussian fixed point: y_t = 2, y_h = 1 + d/2.  At the Ising upper critical
    # dimension d = 4 this gives the classical (Landau) exponents AND — since
    # hyperscaling holds exactly at d_upper — still satisfies Josephson.
    s = ScalingDimensions(2 // 1, 1 + 4 // 2, 4)   # y_h = 3, d = 4
    e = critical_exponents(s)
    @test e == (α=0 // 1, β=1 // 2, γ=1 // 1, δ=3 // 1, ν=1 // 2, η=0 // 1)
    @test exponents_consistent(e; d=4)
end

@testset "inverse map: (ν, η, d) reconstructs the eigenvalues and the rest" begin
    for (yt, yh, dd) in _RG_SETS
        e = critical_exponents(ScalingDimensions(yt, yh, dd))
        # recover the eigenvalues from just ν and η at dimension d …
        s = scaling_dimensions(; ν=e.ν, η=e.η, d=dd)
        @test s.y_t == yt
        @test s.y_h == yh
        # … and the FULL exponent set round-trips (δ, β, γ, α all reconstructed
        # from ν, η, d alone — the two-eigenvalue structure)
        @test critical_exponents(s) == e
    end
end

@testset "float eigenvalues propagate (numerical fixed point)" begin
    s = ScalingDimensions(1.0, 1.875, 2.0)
    e = critical_exponents(s)
    @test e.β ≈ 0.125
    @test e.γ ≈ 1.75
    @test e.δ ≈ 15.0
    @test e.η ≈ 0.25
    # laws still hold to floating tolerance
    @test residual(Rushbrooke(); α=e.α, β=e.β, γ=e.γ) ≈ 0 atol = 1e-12
    @test residual(Fisher(); γ=e.γ, ν=e.ν, η=e.η) ≈ 0 atol = 1e-12
end
