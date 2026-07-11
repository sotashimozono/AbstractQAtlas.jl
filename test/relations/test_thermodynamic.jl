# FDT identities: fluctuation route vs INDEPENDENT derivative route.
#
# For closed-form systems (two-level, two free spins) the moments come
# from the exact ensemble while the response comes from a numerical
# derivative of the exact ⟨·⟩(T or h) — two genuinely different
# computations of the same physical response.  Central differences with
# step h have O(h²) truncation error; the tolerances below follow.

using AbstractQAtlas
using AbstractQAtlas: residual, check, solve

# exact two-level system: E ∈ {0, ε}
tl_Z(β, ε) = 1 + exp(-β * ε)
tl_E(β, ε) = ε * exp(-β * ε) / tl_Z(β, ε)          # ⟨E⟩
tl_varE(β, ε) = ε^2 * exp(-β * ε) / tl_Z(β, ε)^2   # ⟨E²⟩ − ⟨E⟩²

@testset "SpecificHeatFDT: β²Var(E) == dE/dT (two-level, N = 1)" begin
    ε = 1.3
    for T in (0.3, 0.7, 1.5, 4.0)
        β = 1 / T
        # response route: C = d⟨E⟩/dT by central difference of the exact ⟨E⟩(T)
        h = 1e-4 * T
        C_deriv = (tl_E(1 / (T + h), ε) - tl_E(1 / (T - h), ε)) / (2h)
        # fluctuation route: the relation's solve = the estimator
        C_fluct = solve(SpecificHeatFDT(), Val(:C); var_E=tl_varE(β, ε), β=β)
        @test isapprox(C_fluct, C_deriv; rtol=1e-6)
        @test check(SpecificHeatFDT(); C=C_deriv, var_E=tl_varE(β, ε), β=β, atol=1e-6)
        # T-form kwarg equals β-form
        @test solve(SpecificHeatFDT(), Val(:C); var_E=tl_varE(β, ε), T=T) ≈ C_fluct
        # var_E back-solve round-trips
        @test solve(SpecificHeatFDT(), Val(:var_E); C=C_fluct, β=β) ≈ tl_varE(β, ε)
    end
end

# two independent Ising spins σᵢ = ±1 in field h: H = −h(σ₁+σ₂), M = σ₁+σ₂
ts_M(β, h) = 2 * tanh(β * h)                       # exact ⟨M⟩
ts_varM(β, h) = 2 / cosh(β * h)^2                  # exact Var(M)

@testset "SusceptibilityFDT: βVar(M)/N == (∂⟨M⟩/∂h)/N (two spins)" begin
    N = 2
    for T in (0.5, 1.0, 2.5), h0 in (0.0, 0.4)
        β = 1 / T
        δ = 1e-5
        # response route: χ_per_site = (1/N) d⟨M⟩/dh at h0, central difference
        χ_deriv = (ts_M(β, h0 + δ) - ts_M(β, h0 - δ)) / (2δ) / N
        # fluctuation route
        χ_fluct = solve(SusceptibilityFDT(), Val(:χ); var_M=ts_varM(β, h0), β=β, N=N)
        @test isapprox(χ_fluct, χ_deriv; rtol=1e-6, atol=1e-10)
        @test check(
            SusceptibilityFDT(); χ=χ_deriv, var_M=ts_varM(β, h0), β=β, N=N, atol=1e-8
        )
    end
end

@testset "LinearResponseFDT: ∂⟨O⟩/∂λ == βVar(O) (diagonal observable)" begin
    # H(λ) = H₀ − λO on a random diagonal spectrum: exact Gibbs moments.
    E0 = [0.0, 0.35, 0.8, 1.1, 2.3]
    O = [-2.0, -1.0, 0.0, 1.0, 2.0]
    gibbs(β, λ) = exp.(-β .* (E0 .- λ .* O)) ./ sum(exp.(-β .* (E0 .- λ .* O)))
    meanO(β, λ) = sum(gibbs(β, λ) .* O)
    varO(β, λ) = sum(gibbs(β, λ) .* O .^ 2) - meanO(β, λ)^2
    for β in (0.4, 1.0, 2.0), λ in (0.0, 0.3)
        δ = 1e-5
        dO = (meanO(β, λ + δ) - meanO(β, λ - δ)) / (2δ)
        @test check(LinearResponseFDT(); dO_dλ=dO, var_O=varO(β, λ), β=β, atol=1e-7)
        @test solve(LinearResponseFDT(), Val(:dO_dλ); var_O=varO(β, λ), β=β) ≈ dO rtol =
            1e-6
    end
end

@testset "SpecificHeatFromEntropy: c = T ∂s/∂T vs the closed-form ∂s/∂T" begin
    # two-level system entropy s(T) in closed form; ∂s/∂T by finite diff,
    # and c_v independently from Var(E)/T² — the two routes must agree.
    ε = 1.1
    for T in (0.4, 1.0, 2.5)
        β = 1 / T
        Z = 1 + exp(-β * ε)
        s(TT) =
            (bb=1 / TT; log(1 + exp(-bb * ε)) + bb * ε * exp(-bb * ε) / (1 + exp(-bb * ε)))
        h = 1e-5 * T
        dS_dT = (s(T + h) - s(T - h)) / (2h)
        varE = ε^2 * exp(-β * ε) / Z^2
        C = solve(SpecificHeatFDT(), Val(:C); var_E=varE, β=β)     # fluctuation route
        @test check(SpecificHeatFromEntropy(); C=C, dS_dT=dS_dT, T=T, atol=1e-6)
        @test solve(SpecificHeatFromEntropy(), Val(:C); dS_dT=dS_dT, T=T) ≈ C rtol = 1e-5
    end
end

@testset "HeatCapacityDifference: Mayer relation c_p − c_v = T v α²/κ_T" begin
    # exact by construction: pick α, κ_T, and set the difference accordingly
    T, v, α, κT = 2.0, 1.0, 0.5, 1.3
    diff = T * v * α^2 / κT
    Cv = 1.5
    Cp = Cv + diff
    @test check(HeatCapacityDifference(); Cp=Cp, Cv=Cv, T=T, v=v, α=α, κT=κT, atol=1e-13)
    # c_p ≥ c_v always (κ_T > 0)
    @test Cp >= Cv
    # solve for c_p given the rest
    @test solve(HeatCapacityDifference(), Val(:Cp); Cv=Cv, T=T, v=v, α=α, κT=κT) ≈ Cp
    @test !check(HeatCapacityDifference(); Cp=Cv, Cv=Cv, T=T, v=v, α=α, κT=κT)  # ignores diff
end

@testset "StructureFactorSusceptibility: χ = β S(q→0) (static classical FDT)" begin
    # two independent spins in a field: χ (per site) and S(q=0) = Var(M)/N
    # both come from the same fluctuation, so χ = β S(0) exactly.
    N = 2
    for T in (0.5, 1.0, 2.5), h in (0.0, 0.3)
        β = 1 / T
        varM = 2 / cosh(β * h)^2                # exact Var(M) for two spins
        Sq0 = varM / N                          # static structure factor at q=0
        χ = solve(SusceptibilityFDT(), Val(:χ); var_M=varM, β=β, N=N)
        @test check(StructureFactorSusceptibility(); χ=χ, Sq0=Sq0, β=β, atol=1e-13)
        @test solve(StructureFactorSusceptibility(), Val(:χ); Sq0=Sq0, β=β) ≈ χ
    end
end

@testset "argument validation" begin
    @test_throws ErrorException residual(SpecificHeatFDT(); C=1.0, var_E=1.0)  # no β, no T
end
