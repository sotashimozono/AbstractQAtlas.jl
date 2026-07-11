# CFT finite-size relations: read the central charge and scaling
# dimensions off synthetic finite-size spectra, the way an MPS/ED
# calculation would.  Independent check = the extracted universal number
# is the SAME across system sizes (the hallmark of a universal amplitude).

using AbstractQAtlas
using AbstractQAtlas: residual, check, solve

@testset "CasimirCentralCharge: e₀(L) = e_∞ − πcv/(6L²) reads off c" begin
    c, v = 1 / 2, 1.0                      # 2D Ising universality
    e∞ = -1.2732
    e0(L) = e∞ - π * c * v / (6 * L^2)     # synthetic finite-size energy density
    for L in (16, 32, 64, 128)
        dE = e0(L) - e∞                    # the finite-size correction
        @test check(CasimirCentralCharge(); dE=dE, c=c, v=v, L=L, atol=1e-13)
        # c recovered identically at every size (universal amplitude)
        @test solve(CasimirCentralCharge(), Val(:c); dE=dE, v=v, L=L) ≈ c
    end
    # a wrong c fails
    L = 32
    @test !check(CasimirCentralCharge(); dE=e0(L) - e∞, c=1.0, v=v, L=L)
end

@testset "FiniteSizeGap: E_x(L) − E₀(L) = 2πvx/L reads off the scaling dimension" begin
    v = 1.3
    for x in (1 / 8, 1.0, 15 / 8)          # Ising σ, ε, and a descendant
        gaps = Dict(L => 2π * v * x / L for L in (16, 32, 64))
        for (L, gap) in gaps
            @test check(FiniteSizeGap(); gap=gap, x=x, v=v, L=L, atol=1e-12)
            @test solve(FiniteSizeGap(), Val(:x); gap=gap, v=v, L=L) ≈ x
        end
        # the extracted dimension is L-independent (universal)
        xs = [solve(FiniteSizeGap(), Val(:x); gap=gaps[L], v=v, L=L) for L in (16, 32, 64)]
        @test all(≈(x), xs)
    end
end

@testset "one-call sweep: finite-size data → central charge + dimension" begin
    # a bag of measured finite-size numbers picks up exactly the cft relations
    data = (
        dE=(-π * 0.5 * 1.0 / (6 * 64^2)),
        c=0.5,
        v=1.0,
        L=64,
        gap=2π * 1.0 * 0.125 / 64,
        x=0.125,
    )
    rep = relation_report(data; atol=1e-12, domain=:cft)
    @test length(rep) == 2
    @test all(row -> row.pass, rep)
end
